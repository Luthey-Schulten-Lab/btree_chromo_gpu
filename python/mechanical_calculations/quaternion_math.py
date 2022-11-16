import numpy as np

# function to multiply two quaternions
def qMult(q1, q2):

    # q1 - quaternion 1
    # q2 - quaternion 2

    qm = np.zeros(4, dtype=np.double)

    qm[0] = q1[0]*q2[0]-np.dot(q1[1:], q2[1:])
    qm[1:] = q1[0]*q2[1:]+q2[0]*q1[1:]+np.cross(q1[1:], q2[1:])

    return qm

# function to calculate conjugate of a quaternion
def qConj(q):

    qC = q

    qC[1:] = -qC[1:]

    return qC

# function to calculate the dot product of a quaternion
def qDot(q):

    return np.dot(q)

# function calculate the inverse of a quaternion
def qInv(q):

    return qConj(q)/qDot(q)

# function to calculate the norm of a quaternion
def qNorm(q):

    return np.linalg.norm(q)

# function to convert a vector to a quaternion
def vec_to_quat(v):

    q = np.zeros(4, dtype=np.double)

    q[1:] = v

    return q

# function to get vector part of quaternion
def quat_to_vec(q):

    return q[1:]

# function to rotate a vector using a quaternion
def rot_vec_by_quat(v,q):

    u_q = vec_to_quat(v)

    qi = qInv(q)

    u_q = qMult(u_q,qi)

    u_q = qMult(q,u_q)

    u = quat_to_vec(u_q)

    return u
