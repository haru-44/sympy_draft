def gaussian_reduce(w, a, b):
    r"""
    Returns a reduced solution `(x, z)` to the congruence
    `X^2 - aZ^2 \equiv 0 \ (mod \ b)` so that `x^2 + |a|z^2` is minimal.

    Details
    =======

    Here ``w`` is a solution of the congruence `x^2 \equiv a \ (mod \ b)`

    References
    ==========

    .. [1] Gaussian lattice Reduction [online]. Available:
           https://web.archive.org/web/20201021115213/http://home.ie.cuhk.edu.hk/~wkshum/wordpress/?p=404
    .. [2] Cremona, J. E., Rusin, D. (2003). Efficient Solution of Rational Conics.
           Mathematics of Computation, 72(243), 1417-1441.
           https://doi.org/10.1090/S0025-5718-02-01480-1

    """
    def _dot(u, v):
        r"""
        Returns a special dot product of the vectors `u = (u_{1}, u_{2})` and
        `v = (v_{1}, v_{2})` which is defined in order to reduce solution of
        the congruence equation `X^2 - aZ^2 \equiv 0 \ (mod \ b)`.
        """
        u_1, u_2 = u
        v_1, v_2 = v
        return (w*u_1 + b*u_2)*(w*v_1 + b*v_2) + abs(a)*u_1*v_1

    u = (0, 1)
    v = (1 if b*w >= 0 else -1, 0)
    if _dot(u, u) < _dot(v, v):
        u, v = v, u
    while _dot(u, u) > _dot(v, v):
        k = _dot(u, v) // _dot(v, v)
        u, v = v, (u[0] - k*v[0], u[1] - k*v[1])
    c = (v[0] - u[0], v[1] - u[1])
    if 2*_dot(u, v) < _dot(u, u) or _dot(u, u) < _dot(c, c):
        c = u
    return c[0]*w + b*c[1], c[0]
