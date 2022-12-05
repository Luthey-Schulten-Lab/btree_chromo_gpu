# =============================================
# Author: Benjamin R. Gilbert
# Email: brg4@illinois.edu
# =============================================

import numpy as np

# function to multiply two quaternions
def qMult(q1, q2):

    # q1 - quaternion 1
    # q2 - quaternion 2

    qm = np.zeros(4, dtype=np.double)

    qm[0] = q1[0]*q2[0]-np.dot(q1[1:], q2[1:])
    qm[1:] = q1[0]*q2[1:]+q2[0]*q1[1:]+np.cross(q1[1:], q2[1:])

    return qm

# function to generate a set of orientation quaternions about a circular ring in a plane
def genQset(N, ntwists, dt, u):

    # N - number of points around ring
    # ntwists - number of twists
    # dt - angle offset
    # u - axis of rotation

    q = np.zeros((N, 4), dtype=np.double)

    t = ntwists*(2*np.pi)*np.arange(N, dtype=np.double)/N+dt

    q[:, 1:] = np.copy(u)

    q[:, 0] = np.cos(t/2.0)
    q[:, 1:] = q[:, 1:]*np.sin(t/2.0)[:, np.newaxis]

    return q

# function to generate a circular ring of monomers in a plane
def initialize_ring(r, R, compaction, ntwists, a, b):

    # r - radius of monomers
    # R - circle radius
    # compaction - factor for compaction of monomers in circle from equilibrium distance
    # ntwists - number of twists around ring

    # number of monomers
    n = int(compaction*np.floor(2*np.pi*R/(2*r)))

    theta = np.arange(n)*(2*np.pi/n)

    if ntwists <= -1:
        ntwists = np.rint(n*(10.8-10.0)/10.8)

    positions = np.zeros((n, 3), dtype=np.double)
    positions[:, 0] = a*R*np.cos(theta)
    positions[:, 1] = b*R*np.sin(theta)

    I = np.eye(3, dtype=np.double)
    ex = I[0, :]
    ey = I[1, :]
    ez = I[2, :]

    q1 = genQset(n, 1, np.pi/2, ez[np.newaxis, :])

    q1C = -np.copy(q1)
    q1C[:, 0] = -1*q1C[:, 0]

    u = np.zeros((n, 4), dtype=np.double)

    u[:, 1:] = ex[np.newaxis, :]

    for i in range(n):

        u[i, :] = qMult(u[i, :], q1C[i, :])
        u[i, :] = qMult(q1[i, :], u[i, :])

    u = u[:, 1:]

    q2 = genQset(n, ntwists, 0.0, u)

    oris = np.zeros((n, 4), dtype=np.double)

    for i in range(n):

        oris[i, :] = qMult(q2[i, :], q1[i, :])

    return n, positions, oris

# function to generate a circular ring of monomers in a plane
def initialize_toroid(r, R, compaction, ntwists, a, b, toroid_R, toroid_ntwists):

    # r - radius of monomers
    # R - circle radius
    # compaction - factor for compaction of monomers in circle from equilibrium distance
    # ntwists - number of twists around ring

    # number of monomers
    n = int(compaction*np.floor(2*np.pi*R/(2*r)))

    theta = np.arange(n)*(2*np.pi/n)

    if ntwists <= -1:
        ntwists = np.rint(n*(10.8-10.0)/10.8)

    positions = np.zeros((n, 3), dtype=np.double)
    positions[:, 0] = (a*R+toroid_R*np.cos(toroid_ntwists*theta))*np.cos(theta)
    positions[:, 1] = (b*R+toroid_R*np.cos(toroid_ntwists*theta))*np.sin(theta)
    positions[:, 2] = toroid_R*np.sin(toroid_ntwists*theta)

    I = np.eye(3, dtype=np.double)
    ex = I[0, :]
    ey = I[1, :]
    ez = I[2, :]

    q1 = genQset(n, 1, np.pi/2, ez[np.newaxis, :])

    q1C = -np.copy(q1)
    q1C[:, 0] = -1*q1C[:, 0]

    u = np.zeros((n, 4), dtype=np.double)

    u[:, 1:] = ex[np.newaxis, :]

    for i in range(n):

        u[i, :] = qMult(u[i, :], q1C[i, :])
        u[i, :] = qMult(q1[i, :], u[i, :])

    u = u[:, 1:]

    q2 = genQset(n, ntwists, 0.0, u)

    oris = np.zeros((n, 4), dtype=np.double)

    for i in range(n):

        oris[i, :] = qMult(q2[i, :], q1[i, :])

    return n, positions, oris

# function to generate a circular ring of monomers in a plane
def initialize_trefoil(r, R, compaction, ntwists):

    # r - radius of monomers
    # R - circle radius
    # compaction - factor for compaction of monomers in circle from equilibrium distance
    # ntwists - number of twists around ring

    # number of monomers
    n = int(compaction*np.floor(2*np.pi*R/(2*r)))

    theta = np.arange(n)*(2*np.pi/n)

    if ntwists <= -1:
        ntwists = np.rint(n*(10.8-10.0)/10.8)

    positions = np.zeros((n, 3), dtype=np.double)
    positions[:, 0] = R*(np.cos(theta)+2*np.cos(2*theta))
    positions[:, 1] = R*(np.sin(theta)-2*np.sin(2*theta))
    positions[:, 2] = -R*np.sin(3*theta)

    I = np.eye(3, dtype=np.double)
    ex = I[0, :]
    ey = I[1, :]
    ez = I[2, :]

    q1 = genQset(n, 1, np.pi/2, ez[np.newaxis, :])

    q1C = -np.copy(q1)
    q1C[:, 0] = -1*q1C[:, 0]

    u = np.zeros((n, 4), dtype=np.double)

    u[:, 1:] = ex[np.newaxis, :]

    for i in range(n):

        u[i, :] = qMult(u[i, :], q1C[i, :])
        u[i, :] = qMult(q1[i, :], u[i, :])

    u = u[:, 1:]

    q2 = genQset(n, ntwists, 0.0, u)

    oris = np.zeros((n, 4), dtype=np.double)

    for i in range(n):

        oris[i, :] = qMult(q2[i, :], q1[i, :])

    return n, positions, oris

def initialize_linked_rings(r, c_s, compaction, ntwists, a, b):

    # r - radius of monomers
    # c_s - scaling for circle radius
    # compaction - factor for compaction of monomers in circle from equilibrium distance
    # ntwists - number of twists around ring

    # create the first catenane
    n_1, positions_1, oris_1 = initialize_ring(r, c_s, compaction, ntwists, a, b)


    # create second catenane of equal size to the first
    n_2 = n_1
    positions_2 = np.copy(positions_1)
    oris_2 = np.copy(oris_1)

    # rotate the second catenane
    theta_cat = np.pi/2.0

    R_x = np.array([[1,0,0],
                    [0,np.cos(theta_cat),-np.sin(theta_cat)],
                    [0,np.sin(theta_cat),np.cos(theta_cat)]])

    positions_2 = np.matmul(R_x,positions_2.T).T

    # translate the two catenanes such that they are linked
    dx = np.zeros((3),dtype=np.single)
    dx[0] = c_s*r*a*(4.0/5.0)
    dx[0] = c_s*r*a*0.5
    positions_1 += dx
    positions_2 -= dx

    # create a arrays with the combined information for the two catenanes
    n = n_1 + n_2
    positions = np.vstack([positions_1,positions_2])
    oris = np.vstack([oris_1,oris_2])


    return n, positions, oris
