# Python version

MATCH = 3
MISMATCH = -3
GAP = -2
GAP_OPEN = -5
GAP_EXT = -1
NEG_INF = -10**12

def _sub_score(a: str, b: str) -> int:
    return MATCH if a == b else MISMATCH

def global_align_linear(q: str, t: str):
    """Needleman–Wunsch, linear gap = -2. Tie-break: diag > up > left."""
    n, m = len(q), len(t)
    dp = [[0]*(m+1) for _ in range(n+1)]
    ptr = [[-1]*(m+1) for _ in range(n+1)]  # 0=diag,1=up,2=left

    for i in range(1, n+1):
        dp[i][0] = i*GAP
        ptr[i][0] = 1
    for j in range(1, m+1):
        dp[0][j] = j*GAP
        ptr[0][j] = 2

    for i in range(1, n+1):
        qi = q[i-1]
        for j in range(1, m+1):
            tj = t[j-1]
            diag = dp[i-1][j-1] + _sub_score(qi, tj)
            up   = dp[i-1][j] + GAP
            left = dp[i][j-1] + GAP
            best, d = diag, 0
            if up > best:
                best, d = up, 1
            if left > best:
                best, d = left, 2
            dp[i][j] = best
            ptr[i][j] = d

    i, j = n, m
    aq, at = [], []
    while i > 0 or j > 0:
        d = ptr[i][j] if i >= 0 and j >= 0 else -1
        if d == 0:
            aq.append(q[i-1]); at.append(t[j-1]); i -= 1; j -= 1
        elif d == 1:
            aq.append(q[i-1]); at.append('-'); i -= 1
        elif d == 2:
            aq.append('-'); at.append(t[j-1]); j -= 1
        else:
            break
    aq_str = ''.join(reversed(aq))
    at_str = ''.join(reversed(at))
    return dp[n][m], aq_str, at_str

# stubs for the rest (we’ll fill these next)
def local_align_linear(q: str, t: str):
    n, m = len(q), len(t)
    H = [[0] * (m + 1) for _ in range(n + 1)]
    ptr = [[-1] * (m + 1) for _ in range(n + 1)]  # -1 = stop; 0=diag; 1=up; 2=left

    best_score = 0
    best_i = 0
    best_j = 0

    for i in range(1, n + 1):
        qi = q[i - 1]
        Hi_1 = H[i - 1]
        Hi = H[i]
        ptri = ptr[i]
        for j in range(1, m + 1):
            tj = t[j - 1]
            diag = Hi_1[j - 1] + (MATCH if qi == tj else MISMATCH)
            up   = Hi_1[j] + GAP      # gap in target
            left = Hi[j - 1] + GAP    # gap in query

            # include 0 for local reset; deterministic with strict '>'
            s = 0
            d = -1
            if diag > s:
                s, d = diag, 0
            if up > s:
                s, d = up, 1
            if left > s:
                s, d = left, 2

            Hi[j] = s
            ptri[j] = d

            # keep the first occurrence (row-major) of the global max
            if s > best_score:
                best_score = s
                best_i, best_j = i, j

    # Traceback from best cell until score hits 0 or ptr == -1
    i, j = best_i, best_j
    aq, at = [], []
    while i > 0 and j > 0:
        if H[i][j] == 0:
            break
        d = ptr[i][j]
        if d == 0:          # diag
            aq.append(q[i - 1]); at.append(t[j - 1]); i -= 1; j -= 1
        elif d == 1:        # up (gap in target)
            aq.append(q[i - 1]); at.append('-'); i -= 1
        elif d == 2:        # left (gap in query)
            aq.append('-'); at.append(t[j - 1]); j -= 1
        else:
            break

    aq_str = ''.join(reversed(aq))
    at_str = ''.join(reversed(at))
    return best_score, aq_str, at_str


def global_align_affine(q: str, t: str):
    n, m = len(q), len(t)

    # Matrices: M (match/mismatch), X (gap in target, i.e., up), Y (gap in query, i.e., left)
    M = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
    X = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
    Y = [[NEG_INF] * (m + 1) for _ in range(n + 1)]

    # Backpointers
    # ptrM: 0=from M diag, 1=from X diag, 2=from Y diag
    # ptrX: 0=from M (open, up), 1=from X (extend, up)
    # ptrY: 0=from M (open, left), 1=from Y (extend, left)
    ptrM = [[-1] * (m + 1) for _ in range(n + 1)]
    ptrX = [[-1] * (m + 1) for _ in range(n + 1)]
    ptrY = [[-1] * (m + 1) for _ in range(n + 1)]

    # Init
    M[0][0] = 0
    for i in range(1, n + 1):
        open_val  = M[i - 1][0] + GAP_OPEN
        extend_val = X[i - 1][0] + GAP_EXT
        if open_val >= extend_val:    # prefer opening on tie
            X[i][0] = open_val;  ptrX[i][0] = 0
        else:
            X[i][0] = extend_val; ptrX[i][0] = 1
    for j in range(1, m + 1):
        open_val  = M[0][j - 1] + GAP_OPEN
        extend_val = Y[0][j - 1] + GAP_EXT
        if open_val >= extend_val:
            Y[0][j] = open_val;  ptrY[0][j] = 0
        else:
            Y[0][j] = extend_val; ptrY[0][j] = 1

    # Fill
    for i in range(1, n + 1):
        qi = q[i - 1]
        for j in range(1, m + 1):
            tj = t[j - 1]
            s = MATCH if qi == tj else MISMATCH

            # M
            m_from_M = M[i - 1][j - 1] + s
            m_from_X = X[i - 1][j - 1] + s
            m_from_Y = Y[i - 1][j - 1] + s
            bestM, srcM = m_from_M, 0
            if m_from_X > bestM:
                bestM, srcM = m_from_X, 1
            if m_from_Y > bestM:
                bestM, srcM = m_from_Y, 2
            M[i][j] = bestM; ptrM[i][j] = srcM

            # X (up)
            x_from_M = M[i - 1][j] + GAP_OPEN
            x_from_X = X[i - 1][j] + GAP_EXT
            if x_from_M >= x_from_X:   # prefer opening on tie
                X[i][j] = x_from_M; ptrX[i][j] = 0
            else:
                X[i][j] = x_from_X; ptrX[i][j] = 1

            # Y (left)
            y_from_M = M[i][j - 1] + GAP_OPEN
            y_from_Y = Y[i][j - 1] + GAP_EXT
            if y_from_M >= y_from_Y:
                Y[i][j] = y_from_M; ptrY[i][j] = 0
            else:
                Y[i][j] = y_from_Y; ptrY[i][j] = 1

    # Choose terminal state at (n,m) with diag > up > left
    state = 0
    best = M[n][m]
    if X[n][m] > best:
        state, best = 1, X[n][m]
    if Y[n][m] > best:
        state, best = 2, Y[n][m]

    # Traceback
    i, j = n, m
    aq, at = [], []
    while i > 0 or j > 0:
        if i == 0 and j > 0:
            state = 2  # must move left
        elif j == 0 and i > 0:
            state = 1  # must move up

        if state == 0:      # M: diag
            src = ptrM[i][j]
            aq.append(q[i - 1]); at.append(t[j - 1])
            i -= 1; j -= 1
            state = src
        elif state == 1:    # X: up (gap in target)
            src = ptrX[i][j]
            aq.append(q[i - 1]); at.append('-')
            i -= 1
            state = 0 if src == 0 else 1
        else:               # Y: left (gap in query)
            src = ptrY[i][j]
            aq.append('-'); at.append(t[j - 1])
            j -= 1
            state = 0 if src == 0 else 2

    aq_str = ''.join(reversed(aq))
    at_str = ''.join(reversed(at))
    return best, aq_str, at_str


def fitting_align_affine(q: str, t: str):
    n, m = len(q), len(t)

    # States
    M = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
    X = [[NEG_INF] * (m + 1) for _ in range(n + 1)]  # up: gap in target
    Y = [[NEG_INF] * (m + 1) for _ in range(n + 1)]  # left: gap in query

    # Backpointers
    # ptrM: 0=from M diag, 1=from X diag, 2=from Y diag
    # ptrX: 0=from M (open, up), 1=from X (extend, up)
    # ptrY: 0=from M (open, left), 1=from Y (extend, left)
    ptrM = [[-1] * (m + 1) for _ in range(n + 1)]
    ptrX = [[-1] * (m + 1) for _ in range(n + 1)]
    ptrY = [[-1] * (m + 1) for _ in range(n + 1)]

    # Init: free target prefix ⇒ top row can start anywhere at 0 via M
    M[0][0] = 0
    for j in range(1, m + 1):
        M[0][j] = 0  # free to enter at any target position
    # First column: allow leading gaps in target (penalized normally)
    for i in range(1, n + 1):
        open_val = M[i - 1][0] + GAP_OPEN
        extend_val = X[i - 1][0] + GAP_EXT
        if open_val >= extend_val:    # prefer opening on tie
            X[i][0] = open_val;  ptrX[i][0] = 0
        else:
            X[i][0] = extend_val; ptrX[i][0] = 1

    # Fill
    for i in range(1, n + 1):
        qi = q[i - 1]
        for j in range(1, m + 1):
            tj = t[j - 1]
            s = MATCH if qi == tj else MISMATCH

            # M (diag from any state)
            m_from_M = M[i - 1][j - 1] + s
            m_from_X = X[i - 1][j - 1] + s
            m_from_Y = Y[i - 1][j - 1] + s
            bestM, srcM = m_from_M, 0
            if m_from_X > bestM:
                bestM, srcM = m_from_X, 1
            if m_from_Y > bestM:
                bestM, srcM = m_from_Y, 2
            M[i][j] = bestM; ptrM[i][j] = srcM

            # X (up)
            x_from_M = M[i - 1][j] + GAP_OPEN
            x_from_X = X[i - 1][j] + GAP_EXT
            if x_from_M >= x_from_X:
                X[i][j] = x_from_M; ptrX[i][j] = 0
            else:
                X[i][j] = x_from_X; ptrX[i][j] = 1

            # Y (left)
            y_from_M = M[i][j - 1] + GAP_OPEN
            y_from_Y = Y[i][j - 1] + GAP_EXT
            if y_from_M >= y_from_Y:
                Y[i][j] = y_from_M; ptrY[i][j] = 0
            else:
                Y[i][j] = y_from_Y; ptrY[i][j] = 1

    # Termination: best across last row (i=n), leftmost j on tie; apply diag>up>left per j
    best_score = NEG_INF
    best_state = 0
    best_j = 0
    i = n
    for j in range(0, m + 1):
        sM, sX, sY = M[i][j], X[i][j], Y[i][j]
        state_j, s_best_j = 0, sM
        if sX > s_best_j:
            state_j, s_best_j = 1, sX
        if sY > s_best_j:
            state_j, s_best_j = 2, sY
        if s_best_j > best_score:
            best_score, best_state, best_j = s_best_j, state_j, j

    # Traceback from (n, best_j, best_state) until query is fully consumed (i==0)
    j = best_j
    state = best_state
    aq, at = [], []
    while i > 0:
        if j == 0 and state == 2:
            state = 1  # cannot move left at j==0; must go up
        if state == 0:      # M
            src = ptrM[i][j]
            aq.append(q[i - 1]); at.append(t[j - 1])
            i -= 1; j -= 1
            state = src
        elif state == 1:    # X: up (gap in target)
            src = ptrX[i][j]
            aq.append(q[i - 1]); at.append('-')
            i -= 1
            state = 0 if src == 0 else 1
        else:               # Y: left (gap in query)
            src = ptrY[i][j]
            aq.append('-'); at.append(t[j - 1])
            j -= 1
            state = 0 if src == 0 else 2
        if j < 0:
            break

    aq_str = ''.join(reversed(aq))
    at_str = ''.join(reversed(at))
    return best_score, aq_str, at_str
