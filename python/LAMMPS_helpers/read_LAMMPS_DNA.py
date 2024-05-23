import numpy as np
import pickle

def read_traj(filename,in_timesteps=None,slices=None):

    traj = dict()

    with open(filename,'r') as f:

        # initial read to determine atom numbers and numbers of timesteps

        timesteps = np.zeros((1),dtype=np.int32)
        Natoms = np.zeros((1),dtype=np.int32)

        first_frame = True
        
        while True:

            line = f.readline()

            if not line:
                break

            if line.split()[-1] == 'TIMESTEP':

                # determine the timestep of the frame
                line = f.readline()
                t = int(line)

                # determine the number of atoms in the frame
                next(f)
                line = f.readline()
                N = int(line)

                if first_frame == True:
                    timesteps[0] = t
                    Natoms[0] = N
                    first_frame = False
                else:
                    timesteps = np.append(timesteps,t)
                    Natoms = np.append(Natoms,N)

                for _ in range(N + 5):
                    f.readline()


    if in_timesteps == 'Last':
        timesteps = np.array(timesteps[-1],dtype=np.int32)
        Natoms = np.array(Natoms[-1],dtype=np.int32)
    elif in_timesteps == 'FirstLast':
        timesteps = np.array((timesteps[0],timesteps[-1]),dtype=np.int32)
        Natoms = np.array((Natoms[0],Natoms[-1]),dtype=np.int32)
    elif in_timesteps == 'All':
        timesteps = np.array(timesteps,dtype=np.int32)
        Natoms = np.array(Natoms,dtype=np.int32)
    elif in_timesteps == 'Slices':
        temp_timesteps = []
        temp_Natoms = []
        for i in range(len(slices)):
            for j in range(timesteps.shape[0]):
                if timesteps[j] == slices[i]:
                    temp_timesteps.append(timesteps[j])
                    temp_Natoms.append(Natoms[j])
        timesteps = np.array(temp_timesteps,dtype=np.int32)
        Natoms = np.array(temp_Natoms,dtype=np.int32)

    with open(filename,'r') as f:
        
        print(timesteps)
        print(Natoms)

        traj['timesteps'] = timesteps

        Nt = timesteps.shape[0]
        Natoms_max = np.max(Natoms)

        for i_t in range(Nt):

            f.seek(0)

            t = str(timesteps[i_t])
            traj[t] = dict()

            skip_row_count = 9*(i_t+1) + np.sum(Natoms[:i_t])

            print('skip_row_count = '+str(skip_row_count))

            # atom ids
            temp_id = np.loadtxt(f,dtype=np.int32,\
                                 skiprows=skip_row_count,\
                                 usecols=(0),\
                                 max_rows=Natoms[i_t])

            f.seek(0)

            # atom types
            temp_type = np.loadtxt(f,dtype=np.int32,\
                                   skiprows=skip_row_count,\
                                   usecols=(1),\
                                   max_rows=Natoms[i_t])

            f.seek(0)

            # atom positions
            temp_x = np.loadtxt(f,dtype=np.float32,\
                                skiprows=skip_row_count,\
                                usecols=(2,3,4),\
                                max_rows=Natoms[i_t])

            f.seek(0)
            
            # atom orientations
            temp_q = np.loadtxt(f,dtype=np.float32,\
                                skiprows=skip_row_count,\
                                usecols=(5,6,7,8),\
                                max_rows=Natoms[i_t])


            # create dictionaries to hold the DNA, ribo, and bdry information
            traj[t]['Natoms'] = Natoms[i_t]
            traj[t]['DNA'] = dict()
            traj[t]['ribo'] = dict()
            traj[t]['bdry'] = dict()

            DNA_idx = np.argwhere(temp_type>2).flatten()
            N_DNA = DNA_idx.shape[0]
            traj[t]['DNA']['N'] = N_DNA
            offset_DNA_id = temp_id[DNA_idx] - 1
            ribo_idx = np.argwhere(temp_type==2).flatten()
            offset_ribo_id = temp_id[ribo_idx] - N_DNA - 1
            N_ribo = ribo_idx.shape[0]
            traj[t]['ribo']['N'] = N_ribo
            bdry_idx = np.argwhere(temp_type==1).flatten()
            N_bdry = bdry_idx.shape[0]
            traj[t]['bdry']['N'] = N_bdry
            offset_bdry_id = temp_id[bdry_idx] - N_DNA - N_ribo - 1

            # get DNA positions
            traj[t]['DNA']['x'] = np.zeros((N_DNA,3),dtype=np.float32)
            for i in range(N_DNA):
                traj[t]['DNA']['x'][offset_DNA_id[i],:] = \
                    temp_x[DNA_idx[i],:]

            # get DNA orientations
            traj[t]['DNA']['q'] = np.zeros((N_DNA,4),dtype=np.float32)
            for i in range(N_DNA):
                traj[t]['DNA']['q'][offset_DNA_id[i],:] = \
                    temp_q[DNA_idx[i],:]

            # get DNA type
            traj[t]['DNA']['type'] = np.zeros((N_DNA),dtype=np.int32)
            for i in range(N_DNA):
                traj[t]['DNA']['type'][offset_DNA_id[i]] = \
                    temp_type[DNA_idx[i]]

            # get ribo positions
            traj[t]['ribo']['x'] = np.zeros((N_ribo,3),dtype=np.float32)
            for i in range(N_ribo):
                traj[t]['ribo']['x'][offset_ribo_id[i],:] = \
                    temp_x[ribo_idx[i],:]

            # get bdry positions
            traj[t]['bdry']['x'] = np.zeros((N_bdry,3),dtype=np.float32)
            for i in range(N_bdry):
                traj[t]['bdry']['x'][offset_bdry_id[i],:] = \
                    temp_x[bdry_idx[i],:]


    return traj


def write_pickle_traj(filename,traj):

    with open(filename,'wb') as f:

        pickle.dump(traj,f,-1)

    return

def read_pickle_traj(filename):

    with open(filename,'rb') as f:

        traj = pickle.load(f)

    return traj
