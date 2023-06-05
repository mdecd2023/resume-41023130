from pyslvs.metaheuristics import de
from pyslvs.metaheuristics import ALGORITHM, AlgorithmType
import pyslvs
from pyslvs.metaheuristics.utility import ObjFunc

fun = pyslvs.metaheuristics.utility.ObjFunc()
print(dir(fun))
fun1 = ObjFunc()

settings = {"1":1, 'max_gen': 10, 'max_time': 10}
myde = de.Differential(fun, settings)
print("de")
'''
cdef class FMatch(ObjFunc):
    """This class is used to verified kinematics of the linkage mechanism.

    A fast matching method that adds mapping angles to variables.
    """
    cdef bint bfgs_mode, shape_only, use_curvature, full_path, ordered
    cdef int target_count, target_len, input_count, l_base
    cdef list vpoints
    cdef long[:] target_nodes, pivots
    cdef double[:, :, :] target
    cdef EStack exprs
    cdef cset[int] slider
    cdef map[Sym, CCoord] joint_pos
    cdef map[Sym, double] param

    def __cinit__(self, dict mech):
        # mech = {
        #     'expression': List[VPoint],
        #     'input': OrderedDict([((b0, d0), [start, end]), ...]),
        #     'placement': {pt: (x, y, r)},
        #     'target': {pt: [(x0, y0), (x1, y1), ...]},
        #     'same': {pt: match_to_pt},
        #     # Bounds of base link length
        #     'upper': float,
        #     'lower': float,
        #     'shape_only': bool,
        #     'use_curvature': bool,
        #     'full_path': bool,
        # }
        placement = mech.get('placement', {})
        if len(placement) == 0:
            raise ValueError("no grounded joint")
        target = mech.get('target', {})
        if len(target) == 0:
            raise ValueError("no target joint")
        check_set = {len(t) for t in target.values()}
        if len(check_set) != 1:
            raise ValueError("target paths should be in the same size")
        self.target_len = check_set.pop()
        # Change the target paths into memory view
        self.target_count = len(target)
        self.target_nodes = zeros(self.target_count, dtype=int)
        self.target = zeros((self.target_count, self.target_len, 2), dtype=f64)
        same = mech.get('same', {})
        self.shape_only = mech.get('shape_only', False)
        self.use_curvature = mech.get('use_curvature', False)
        self.full_path = mech.get('full_path', False)
        if self.use_curvature:
            self.shape_only = False
        cdef int i, j
        cdef double[:, :] path
        for i, j in enumerate(target):
            self.target_nodes[i] = j
            path = array(target[j], dtype=f64)
            if self.shape_only:
                _norm(path, 1)
            if self.use_curvature:
                path = _path_signature(_curvature(path), 100)
            self.target[i, :, :] = path
        # Expressions (must be readonly)
        self.vpoints = list(mech.get('expression', []))
        self.pivots = array([i for i, vp in enumerate(self.vpoints)
                             if (<VPoint>vp).grounded()], dtype=int)
        self.slider = {i for i, vp in enumerate(self.vpoints)
                       if (<VPoint>vp).is_slider()}
        inputs = OrderedDict(mech.get('input', {}))
        self.input_count = len(inputs)
        status = {}
        self.exprs = t_config(self.vpoints, tuple(inputs.keys()), status)
        self.bfgs_mode = not all(status.values())
        if not preprocessing(self.exprs, self.vpoints, {p: 0. for p in inputs},
                             self.joint_pos, map[SwappablePair, double](),
                             self.param):
            raise ValueError("wrong number of input parameters")
        # Bounds
        ub = []
        lb = []
        # Position
        cdef double x, y, r
        for i in self.pivots:
            x, y, r = placement[i]
            ub.append(x + r)
            ub.append(y + r)
            lb.append(x - r)
            lb.append(y - r)
        # Length of links
        link_upper = float(mech.get('upper', 100))
        link_lower = float(mech.get('lower', 0))
        cdef Expr expr
        for expr in self.exprs.stack:
            ub.append(link_upper)
            lb.append(link_lower)
            if expr.func in {PLA, PLAP} and expr.v2.first == A_LABEL:
                # The included angle of the link
                ub.append(2 * M_PI)
                lb.append(0)
            elif expr.func == PLLP:
                ub.append(link_upper)
                lb.append(link_lower)
        # The start of the angle parameters
        self.l_base = len(ub)
        if self.use_curvature and self.full_path:
            # Scale factor
            ub.append(1)
            lb.append(1e-12)
        else:
            # Input nodes
            for start, end in inputs.values():
                ub.append(start / 180 * M_PI)
                lb.append(end / 180 * M_PI)
            # Angle rage (input count * target count)
            ub[self.l_base:] *= self.target_len
            lb[self.l_base:] *= self.target_len
        self.ub = array(ub, dtype=f64)
        self.lb = array(lb, dtype=f64)
        # Swap upper and lower bound if reversed
        for i in range(len(self.ub)):
            if self.ub[i] < self.lb[i]:
                self.ub[i], self.lb[i] = self.lb[i], self.ub[i]

    cpdef bint is_two_kernel(self):
        """Input a generic data (variable array), return the mechanism
        expression.
        """
        return self.bfgs_mode

    cdef double fitness(self, double[:] v) nogil:
        """Return the difference of the path signature.

        + Position of fix joints.
        + Link lengths.
        + Angle corresponding to the target points.
        """
        cdef map[int, vector[CCoord]] target
        cdef map[Sym, double] param = map[Sym, double](self.param)
        cdef map[Sym, CCoord] joint_pos = map[Sym, CCoord](self.joint_pos)
        cdef bint ok
        cdef int i, j, node, vi
        cdef double x, y
        cdef Expr expr
        cdef CCoord c
        cdef pair[Sym, CCoord] jp
        for i in range(self.target_len):
            vi = 0
            for node in self.pivots:
                joint_pos[Sym(P_LABEL, node)] = CCoord(v[vi], v[vi + 1])
                vi += 2
            # Input parameters (length)
            for expr in self.exprs.stack:
                param[expr.v1] = v[vi]
                vi += 1
                if expr.v2.first == I_LABEL:
                    continue
                param[expr.v2] = v[vi]
                vi += 1
            # Input parameters (angles)
            for vi in range(self.input_count):
                param[Sym(I_LABEL, vi)] = v[self.l_base + vi
                                            + i * self.target_count]
            # Solve
            ok, joint_pos = quick_solve(self.exprs.stack, joint_pos, param)
            if not ok:
                return HUGE_VAL
            if self.bfgs_mode:
                with gil:
                    data_dict = {}
                    for jp in joint_pos:
                        data_dict[jp.first.second] = Coord.__new__(Coord,
                                                                   jp.second.x,
                                                                   jp.second.y)
                    # Solve
                    try:
                        solved_bfgs = SolverSystem(self.vpoints, {}, data_dict).solve()
                    except ValueError:
                        return HUGE_VAL
            # Collecting
            for node in self.target_nodes:
                if self.bfgs_mode:
                    with gil:
                        x, y = solved_bfgs[node][0]
                    target[node].push_back(CCoord(x, y))
                else:
                    c = joint_pos[Sym(P_LABEL, node)]
                    target[node].push_back(c)
        # Compare
        cdef double fitness = 0
        cdef double scale
        cdef double[:, :] path1, path2
        for vi, node in enumerate(self.target_nodes):
            if self.use_curvature:
                with gil:
                    path1 = zeros((self.target_len, 2), dtype=f64)
                    path2 = array(self.target[vi], dtype=f64)
                for i in range(self.target_len):
                    c = target[node][i]
                    path1[i, 0] = c.x
                    path1[i, 1] = c.y
                path1 = _slice_nan2d(path1)
                if len(path1) == 0:
                    return HUGE_VAL
                path1 = _path_signature(_curvature(path1), 100)
                scale = 1 / v[len(v) - 1]
                if not self.full_path:
                    scale *= _extr1d(path1[:, 0], 1) / _extr1d(path2[:, 0], 1)
                _mul1d(path2[:, 0], scale)
                with gil:
                    j = argmax(_cross_correlation(path2, path1, 0.1))
                for i in range(len(path2)):
                    path2[i, 0] += j
                for i in range(self.target_len):
                    fitness += path1[i, 0] - path2[i, 0]
            else:
                if self.shape_only:
                    with gil:
                        path1 = zeros((self.target_len, 2), dtype=f64)
                    for i in range(self.target_len):
                        c = target[node][i]
                        path1[i, 0] = c.x
                        path1[i, 1] = c.y
                    _norm(path1, 1)
                    for i in range(self.target_len):
                        target[node][i].x = path1[i, 0]
                        target[node][i].y = path1[i, 1]
                for i in range(self.target_len):
                    c = target[node][i]
                    fitness += distance(c.x, c.y,
                                        self.target[vi, i, 0],
                                        self.target[vi, i, 1])
        return fitness

    cpdef object result(self, double[:] v):
        """Input a generic data (variable array), return the mechanism
        expression.
        """
        cdef map[Sym, double] param = map[Sym, double](self.param)
        cdef map[Sym, CCoord] joint_pos = map[Sym, CCoord](self.joint_pos)
        cdef int vi = 0
        cdef int node
        for node in self.pivots:
            joint_pos[Sym(P_LABEL, node)] = CCoord(v[vi], v[vi + 1])
            vi += 2
        cdef Expr expr
        for expr in self.exprs.stack:
            param[expr.v1] = v[vi]
            vi += 1
            if expr.v2.first == I_LABEL:
                continue
            param[expr.v2] = v[vi]
            vi += 1
        # Input parameters (angles)
        for vi in range(self.input_count):
            param[Sym(I_LABEL, vi)] = v[self.l_base + vi]
        # Solve
        _, joint_pos = quick_solve(self.exprs.stack, joint_pos, param)
        cdef pair[Sym, CCoord] jp
        if self.bfgs_mode:
            data_dict = {}
            for jp in joint_pos:
                data_dict[jp.first.second] = Coord.__new__(Coord,
                                                           jp.second.x,
                                                           jp.second.y)
            # Solve
            solved_bfgs = SolverSystem(self.vpoints, {}, data_dict).solve()
        # Collecting
        cdef CCoord c
        cdef VPoint vp
        cdef double x, y
        for node in range(len(self.vpoints)):
            vp = self.vpoints[node]
            if self.bfgs_mode:
                x, y = solved_bfgs[node][0]
                vp.locate(x, y)
                if vp.is_slider():
                    vp.move(solved_bfgs[node][0], solved_bfgs[node][1])
            else:
                c = joint_pos[Sym(P_LABEL, node)]
                if vp.is_slider():
                    vp.locate(vp.c[0, 0], vp.c[0, 1])
                    vp.move((vp.c[0, 0], vp.c[0, 1]), (c.x, c.y))
                else:
                    vp.locate(c.x, c.y)
        expr_str = "M[" + ", ".join([(<VPoint> vp).expr()
                                     for vp in self.vpoints]) + "]"
        logger.debug(f"Result: {expr_str}")
        return expr_str
'''
