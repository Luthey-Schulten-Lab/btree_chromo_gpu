import numpy as np


class mono_range:

    def __init__(self, wrapped, size, ll, mid_ll, mid_ul, ul):

        self.wrapped = wrapped
        self.size = size
        self.ll = ll
        self.mid_ll = mid_ll
        self.mid_ul = mid_ul
        self.ul = ul


    def get_size(self):

        return self.size

    def get_mono_array(self):

        if self.wrapped == False:

            a = np.arange(self.ll,self.ul+1,1,dtype=np.int32)

        else:
            
            a_ll = np.arange(self.ll,self.mid_ul+1,1,dtype=np.int32)
            a_ul = np.arange(self.mid_ul,self.ul+1,1,dtype=np.int32)

            a = np.concatenate(a_ll,a_ul)

        return a

class single_fork_partition:

    def __init__(self, fork, left_size, right_size):

        self.fork = fork
        self.left_size = left_size
        self.right_size = right_size
        self.N_left = 0
        self.N_right = 0
        self.left_mono_ranges = []
        self.right_mono_ranges = []

    def get_fork(self):

        return self.fork

    def get_left_size(self):

        return self.left_size

    def get_right_size(self):

        return self.right_size

    def add_mono_range(self, d, wrapped, size, ll, mid_ll, mid_ul, ul):

        temp_mono_range = mono_range(wrapped,
                                     size,
                                     ll,
                                     mid_ll,
                                     mid_ul,
                                     ul)

        if d == 'l':

            self.N_left += 1
            self.left_mono_ranges.append(temp_mono_range)

        elif d == 'r':

            self.N_right += 1
            self.right_mono_ranges.append(temp_mono_range)


    def get_full_left_mono_array(self):

        a = np.zeros((self.left_size),dtype=np.int32)

        c = 0
        
        for i in range(self.N_left):

            m_r = self.left_mono_ranges[i]

            a[c:c+m_r.get_size()] = m_r.get_mono_array()

            c += m_r.get_size()
            
        return a
    
    def get_full_right_mono_array(self):
    
        a = np.zeros((self.right_size),dtype=np.int32)

        c = 0
        
        for i in range(self.N_right):

            m_r = self.right_mono_ranges[i]

            a[c:c+m_r.get_size()] = m_r.get_mono_array()

            c += m_r.get_size()
            
        return a


class fork_partitions:

    def __init__(self):

        self.N_forks = 0
        self.idx = 0
        self.forks = []

    def get_N_forks(self):

        return self.N_forks

    def get_fork_list(self):

        fork_list = []

        for fork in self.forks:

            fork_list.append(fork.get_fork())

        return fork_list

    def parse_f_p_file(self, filename):

        with open(filename, 'r') as f:

            s = f.readline()
            self.N_forks = int(s.split('=')[1])
            s = f.readline()
            self.idx = int(s.split('=')[1])

            for i_fork in range(self.N_forks):

                fork_line = f.readline()
                fork_line = fork_line.split(',')

                fork = fork_line[0]
                N_left = int(fork_line[1])
                left_size = int(fork_line[2])
                N_right = int(fork_line[3])
                right_size = int(fork_line[4])

                temp_s_f_p = single_fork_partition(fork,
                                                   left_size,
                                                   right_size)

                for i_m_r in range(N_left + N_right):

                    m_r_line = f.readline()
                    m_r_line = m_r_line.strip()
                    m_r_line = m_r_line.split(',')

                    if m_r_line[0] == 'w':
                        wrapped = True
                    else:
                        wrapped = False

                    size = int(m_r_line[1])
                    ll = int(m_r_line[2]) - self.idx
                    mid_ll = int(m_r_line[3]) - self.idx
                    mid_ul = int(m_r_line[4]) - self.idx
                    ul = int(m_r_line[5]) - self.idx

                    if i_m_r < N_left:
                        d = 'l'
                    else:
                        d = 'r'
                    
                    temp_s_f_p.add_mono_range(d,
                                              wrapped,
                                              size,
                                              ll,
                                              mid_ll,
                                              mid_ul,
                                              ul)

                self.forks.append(temp_s_f_p)

    def get_fork_partition(self,i):

        return self.forks[i]
