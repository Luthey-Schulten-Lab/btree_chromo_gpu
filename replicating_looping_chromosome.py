"""
ReplicatingLoopingChromosome: A replicating chromosome with SMC-mediated loop extrusion.

Contains a ReplicatingChromosome for structural bonds and manages loop bonds.
"""

import random
import numpy as np
from replicating_chromosome import ReplicatingChromosome


class ReplicatingLoopingChromosome:
    """
    A replicating chromosome with SMC-mediated loop extrusion.
    Contains a ReplicatingChromosome for structural bonds and manages loop bonds.
    """
    
    def __init__(self, initial_length, num_smc=50, smc_width=1, 
                 basal_death_prob=0.001, step_prob=1.0, 
                 stall_death_prob=0.001, knockoff=0.0, bypass=1.0,
                 random_seed=None):
        """
        Initialize a replicating chromosome with loop extrusion.
        
        Args:
            initial_length: Initial chromosome length
            num_smc: Number of SMC complexes (loops)
            smc_width: Width of SMC occupancy (occupies pos-width to pos+width)
            basal_death_prob: Probability of loop dissociation per timestep
            step_prob: Probability of taking a step
            stall_death_prob: Probability of dissociation when stalled
            knockoff: Probability of knockoff vs bypass when blocked
            bypass: Probability of bypass vs knockoff when blocked
            random_seed: Random seed for reproducibility
        """
        # Initialize chromosome (coords can be loaded later if needed)
        self.chromo = ReplicatingChromosome(initial_length)
        self.initial_length = initial_length
        
        # Loop extrusion parameters
        self.num_smc = num_smc
        self.smc_width = smc_width
        self.basal_death_prob = basal_death_prob
        self.step_prob = step_prob
        self.stall_death_prob = stall_death_prob
        self.knockoff = knockoff
        self.bypass = bypass
        
        # Random number generator
        if random_seed is not None:
            random.seed(random_seed)
            np.random.seed(random_seed)
        
        # Loop data structures
        # Each loop has two anchors (loops1[i], loops2[i])
        # These are bead indices in the chromosome
        self.loops1 = [0] * num_smc  # Left anchor positions (bead indices)
        self.loops2 = [0] * num_smc  # Right anchor positions (bead indices)
        self.stalled1 = [0] * num_smc  # Is left anchor stalled?
        self.stalled2 = [0] * num_smc  # Is right anchor stalled?
        
        # Occupancy tracking: how many SMCs occupy each bead
        # Size: max number of beads (2 * initial_length when fully replicated)
        max_beads = 2 * initial_length
        self.occupied = [0] * (max_beads + 1)
        
        # Track current max bead index (like N in C++ code)
        # This gets updated during replication, no need to query chromosome
        self.current_max_bead = initial_length - 1
        
        # Cache for performance
        self._beads_cache_invalid = True
        
        # Initialize loops by birthing them
        for i in range(num_smc):
            self.birth(i)
    
    def get_all_beads(self):
        """Get all bead indices in the chromosome."""
        return self.chromo.get_all_beads()
    
    def wrapped_index(self, idx):
        """Handle wraparound for bead indices. Uses current_max_bead (like N in C++)."""
        max_bead = self.current_max_bead
        if idx < 0:
            return max_bead + idx + 1
        if idx > max_bead:
            return idx - (max_bead + 1)
        return idx
    
    def occupy_around(self, pos):
        """Mark beads around position pos as occupied."""
        for j in range(-self.smc_width, self.smc_width + 1):
            wrapped_pos = self.wrapped_index(pos + j)
            if wrapped_pos < len(self.occupied):
                self.occupied[wrapped_pos] += 1
    
    def release_around(self, pos):
        """Release occupancy around position pos."""
        for j in range(-self.smc_width, self.smc_width + 1):
            wrapped_pos = self.wrapped_index(pos + j)
            if wrapped_pos < len(self.occupied):
                self.occupied[wrapped_pos] = max(0, self.occupied[wrapped_pos] - 1)
    
    def neighbors(self, pos):
        """Get occupancy values for neighbors around position pos. Optimized."""
        # Pre-allocate list
        width_range = 2 * self.smc_width + 1
        result = [0] * width_range
        for idx, j in enumerate(range(-self.smc_width, self.smc_width + 1)):
            wrapped_pos = self.wrapped_index(pos + j)
            if wrapped_pos < len(self.occupied):
                result[idx] = self.occupied[wrapped_pos]
        return result
    
    def is_empty_around(self, pos):
        """Check if area around position pos is empty. Optimized."""
        # Direct check without creating list
        for j in range(-self.smc_width, self.smc_width + 1):
            wrapped_pos = self.wrapped_index(pos + j)
            if wrapped_pos < len(self.occupied) and self.occupied[wrapped_pos] > 0:
                return False
        return True
    
    def _is_replication_in_progress(self):
        """Check if replication has started but not finished."""
        return (self.chromo.replicated_left > 0 or self.chromo.replicated_right > 0) and not self.chromo.is_fully_replicated()
    
    def no_forks_around(self, pos):
        """Check if replication forks are not near position pos. Optimized."""
        # If replication hasn't started or is complete, no forks to block
        if not self._is_replication_in_progress():
            return True
        
        left_fork = self.chromo.get_left_fork_position()
        right_fork = self.chromo.get_right_fork_position()
        
        for j in range(-self.smc_width, self.smc_width + 1):
            wrapped_pos = self.wrapped_index(pos + j)
            if wrapped_pos == left_fork or wrapped_pos == right_fork:
                return False
        return True
    
    def _get_cached_beads(self):
        """Get beads and cache them for performance."""
        if not hasattr(self, '_cached_beads') or self._beads_cache_invalid:
            self._cached_beads = self.get_all_beads()
            self._beads_cache_invalid = False
        return self._cached_beads
    
    def get_region_for_bead(self, bead_idx):
        """
        Determine which region a bead belongs to.
        Returns: 'left_daughter' or 'right_daughter'
        """
        if bead_idx < self.initial_length:
            return 'left_daughter'
        else:
            return 'right_daughter'
    
    def get_next_position(self, pos, direction):
        """
        Get next position based on chromosome topology.
        Uses current_max_bead instead of querying chromosome.
        direction: -1 for left anchor (moves left), +1 for right anchor (moves right)
        """
        if pos < self.initial_length:
            # Left daughter: circular, indices 0 to initial_length-1
            region_size = self.initial_length
            wrapped = (pos + direction) % region_size
            if wrapped < 0:
                wrapped += region_size
            return wrapped
        else:
            # Right daughter: linear during replication, circular when complete
            right_daughter_min = self.initial_length
            right_daughter_max = self.current_max_bead
            region_size = right_daughter_max - right_daughter_min + 1
            
            zero_based_pos = pos - right_daughter_min
            new_zero_based = zero_based_pos + direction
            
            if new_zero_based < 0:
                # Wrap around if fully replicated, otherwise can't go before start
                if self.chromo.is_fully_replicated():
                    # Wrap around: convert negative to positive modulo
                    wrapped_zero_based = new_zero_based % region_size
                    return right_daughter_min + wrapped_zero_based
                else:
                    return pos  # Can't go before start if not fully replicated
            if new_zero_based >= region_size:
                # Wrap around if fully replicated
                if self.chromo.is_fully_replicated():
                    return right_daughter_min + (new_zero_based % region_size)
                else:
                    return pos  # Can't go after end if not fully replicated
            return right_daughter_min + new_zero_based
    
    def birth(self, i):
        """
        Create a new loop at random positions.
        Loops are created between two nearby beads.
        Uses current_max_bead instead of querying chromosome (like C++ code).
        """
        max_attempts = 1000
        spacing = 2 * self.smc_width + 1
        
        for attempt in range(max_attempts):
            # Pick a random position from 0 to current_max_bead (like C++: 1 to N)
            start_pos = random.randint(0, self.current_max_bead)
            
            if start_pos < self.initial_length:
                # Left daughter: circular topology
                idx1 = start_pos
                # Find idx2 that's spacing away, wrapping if needed
                idx2 = (idx1 + spacing) % self.initial_length
            else:
                # Right daughter: linear during replication
                idx1 = start_pos
                idx2 = idx1 + spacing
                # Check if idx2 is within valid range
                if idx2 > self.current_max_bead:
                    continue
            
            # Check if positions are valid
            if not (self.is_empty_around(idx1) and 
                    self.is_empty_around(idx2) and
                    self.no_forks_around(idx1) and
                    self.no_forks_around(idx2)):
                continue
            
            # Accept this position
            self.loops1[i] = idx1
            self.loops2[i] = idx2
            self.occupy_around(idx1)
            self.occupy_around(idx2)
            self.stalled1[i] = 0
            self.stalled2[i] = 0
            return
        
        # If we couldn't find a position after many attempts, leave as 0,0
        print(f"Warning: Could not find valid position for loop {i} after {max_attempts} attempts")
    
    def death(self):
        """Handle loop dissociation based on death probabilities."""
        for i in range(self.num_smc):
            if self.loops1[i] == 0 and self.loops2[i] == 0:
                continue  # Skip inactive loops
            
            # Calculate death probability based on stall state
            falloff1 = self.stall_death_prob if self.stalled1[i] else self.basal_death_prob
            falloff2 = self.stall_death_prob if self.stalled2[i] else self.basal_death_prob
            falloff_rate = max(falloff1, falloff2)
            
            if random.random() < falloff_rate:
                # Dissociate loop
                self.release_around(self.loops1[i])
                self.release_around(self.loops2[i])
                self.stalled1[i] = 0
                self.stalled2[i] = 0
                self.birth(i)  # Rebirth at new location
    
    def get_next_vacancy(self, pos, direction, max_search=15):
        """Find the next empty position in the given direction."""
        candidate = pos
        for _ in range(max_search):
            candidate = self.get_next_position(candidate, direction)
            if self.is_empty_around(candidate) and self.no_forks_around(candidate):
                return candidate
        return pos  # No vacancy found
    
    def step(self, tag_extruded=False):
        """
        Move loops by extruding (left anchor moves left, right anchor moves right).
        Optimized version.
        tag_extruded: whether to track extruded positions (not implemented yet)
        """
        # Cache fork positions for this step (they don't change during a step)
        # Only get fork positions if replication is in progress
        if self._is_replication_in_progress():
            left_fork = self.chromo.get_left_fork_position()
            right_fork = self.chromo.get_right_fork_position()
        else:
            left_fork = None
            right_fork = None
        
        for i in range(self.num_smc):
            if self.loops1[i] == 0 and self.loops2[i] == 0:
                continue  # Skip inactive loops
            
            cur1 = self.loops1[i]  # Left anchor
            cur2 = self.loops2[i]  # Right anchor
            
            # Reset stall flags
            self.stalled1[i] = 0
            self.stalled2[i] = 0
            
            # Move left anchor (direction -1)
            if random.random() <= self.step_prob:
                new_cur1 = self.get_next_position(cur1, -1)
                # Check only the leftmost neighbor (the new position we're moving into)
                # This matches C++ behavior: neighbors(newCur1).front() == 0
                wrapped_pos = self.wrapped_index(new_cur1 - self.smc_width)
                leftmost_neighbor_empty = (wrapped_pos >= len(self.occupied) or self.occupied[wrapped_pos] == 0)
                
                # Check if forks are nearby (forks block SMCs only during active replication)
                no_forks = self._no_forks_at_pos(new_cur1, left_fork, right_fork)
                
                if leftmost_neighbor_empty and no_forks:
                    # Take normal step
                    self.occupy_around(new_cur1)
                    self.release_around(cur1)
                    self.loops1[i] = new_cur1
                else:
                    # Blocked - try bypass or stall
                    # If blocked by fork, always stall
                    if not no_forks:
                        self.stalled1[i] = 1
                    else:
                        rate_sum = self.knockoff + self.bypass
                        if random.random() <= rate_sum:
                            if new_cur1 < len(self.occupied) and self.occupied[new_cur1] > 2:
                                # Stall
                                self.stalled1[i] = 1
                            elif random.random() <= self.knockoff / rate_sum:
                                # Stall
                                self.stalled1[i] = 1
                            else:
                                # Bypass
                                new_cur1 = self.get_next_vacancy(cur1, -1)
                                self.occupy_around(new_cur1)
                                self.release_around(cur1)
                                self.loops1[i] = new_cur1
            
            # Move right anchor (direction +1)
            if random.random() <= self.step_prob:
                new_cur2 = self.get_next_position(cur2, 1)
                # Check only the rightmost neighbor (the new position we're moving into)
                # This matches C++ behavior: neighbors(newCur2).back() == 0
                wrapped_pos = self.wrapped_index(new_cur2 + self.smc_width)
                rightmost_neighbor_empty = (wrapped_pos >= len(self.occupied) or self.occupied[wrapped_pos] == 0)
                
                # Check if forks are nearby (forks block SMCs only during active replication)
                no_forks = self._no_forks_at_pos(new_cur2, left_fork, right_fork)
                
                if rightmost_neighbor_empty and no_forks:
                    # Take normal step
                    self.occupy_around(new_cur2)
                    self.release_around(cur2)
                    self.loops2[i] = new_cur2
                else:
                    # Blocked - try bypass or stall
                    # If blocked by fork, always stall
                    if not no_forks:
                        self.stalled2[i] = 1
                    else:
                        rate_sum = self.knockoff + self.bypass
                        if random.random() <= rate_sum:
                            if new_cur2 < len(self.occupied) and self.occupied[new_cur2] > 2:
                                # Stall
                                self.stalled2[i] = 1
                            elif random.random() <= self.knockoff / rate_sum:
                                # Stall
                                self.stalled2[i] = 1
                            else:
                                # Bypass
                                new_cur2 = self.get_next_vacancy(cur2, 1)
                                self.occupy_around(new_cur2)
                                self.release_around(cur2)
                                self.loops2[i] = new_cur2
    
    def _no_forks_at_pos(self, pos, left_fork, right_fork):
        """Check if position is not near either fork. Uses _is_pos_near_fork."""
        # If replication hasn't started or is complete, no forks to block
        if not self._is_replication_in_progress() or left_fork is None or right_fork is None:
            return True
        return not (self._is_pos_near_fork(pos, left_fork) or self._is_pos_near_fork(pos, right_fork))
    
    def get_loop_bonds(self):
        """
        Get all loop bonds as a list of tuples (bead1, bead2).
        These represent the connections between loop anchors.
        """
        loop_bonds = []
        for i in range(self.num_smc):
            if self.loops1[i] != 0 or self.loops2[i] != 0:  # Active loop
                bond = tuple(sorted([self.loops1[i], self.loops2[i]]))
                loop_bonds.append(bond)
        return loop_bonds
    
    def get_all_bonds(self):
        """
        Get all bonds: structural bonds from chromosome + loop bonds.
        Returns a list of tuples (bead1, bead2) where bead1 < bead2.
        """
        all_bonds = set(self.chromo.get_bonds())  # Structural bonds
        loop_bonds = self.get_loop_bonds()
        all_bonds.update(loop_bonds)  # Add loop bonds
        return sorted(list(all_bonds))
    
    def write_bonds_to_file(self, filename, k_stretch=1000, l0=0.472, k_angle=3, theta0=0, vertex_offset=0):
        """
        Write bonds to file in dna_bonds format.
        Structural bonds come first, then SMC bonds (annotated with ;SMC).
        
        Args:
            filename: Path to output bonds file
            k_stretch: Stretch constant (default: 1000)
            l0: Equilibrium bond length (default: 0.472)
            k_angle: Angle constant (default: 3)
            theta0: Equilibrium angle (default: 0)
            vertex_offset: Offset to add to all DNA bead indices (default: 0)
                          This accounts for membrane vertices that come before DNA beads
        """
        structural_bonds = self.chromo.get_bonds()
        smc_bonds = self.get_loop_bonds()
        
        with open(filename, 'w') as f:
            # Write header
            f.write("; DNA bonds file (generated by dts_generate)\n")
            f.write("; Global bond parameters:\n")
            f.write(f"; k_stretch = {k_stretch}\n")
            f.write(f"; l0 = {l0}\n")
            f.write(f"; k_angle = {k_angle}\n")
            f.write(f"; theta0 = {theta0}\n")
            f.write("; Number of chains: 2\n")
            f.write("; Format: v1 v2 [k_stretch] [l0] (if k_stretch and l0 are in header, only v1 v2 needed)\n")
            f.write("\n")
            
            # Write structural bonds (with offset)
            for bond in structural_bonds:
                bead1, bead2 = bond
                f.write(f"{bead1 + vertex_offset} {bead2 + vertex_offset}\n")
            
            # Write SMC bonds (annotated with ;SMC, with offset)
            for bond in smc_bonds:
                bead1, bead2 = bond
                f.write(f"{bead1 + vertex_offset} {bead2 + vertex_offset} ;SMC\n")
    
    def load_coords_from_file(self, filename, vertex_offset=0):
        """Load coordinates from dna.q file (delegates to chromosome)."""
        self.chromo.load_coords_from_file(filename, vertex_offset=vertex_offset)
    
    def write_coords_to_file(self, filename, vertex_offset=0):
        """Write coordinates to dna.q file (delegates to chromosome)."""
        self.chromo.write_coords_to_file(filename, vertex_offset=vertex_offset)
    
    def _is_pos_near_fork(self, smc_pos, fork_pos):
        """Check if SMC position is within smc_width of fork position (considering topology)."""
        # If replication is complete, forks don't exist, so SMC can't be near a fork
        if self.chromo.is_fully_replicated():
            return False
        
        # After replication starts, forks and SMCs must be in the same chromosome region
        # Left daughter: indices 0 to initial_length-1
        # Right daughter: indices >= initial_length
        smc_region = 'left_daughter' if smc_pos < self.initial_length else 'right_daughter'
        fork_region = 'left_daughter' if fork_pos < self.initial_length else 'right_daughter'
        
        # If SMC and fork are in different regions, they can't be near each other
        if smc_region != fork_region:
            return False
        
        # Check all positions within smc_width of the fork
        for j in range(-self.smc_width, self.smc_width + 1):
            fork_neighbor = self.wrapped_index(fork_pos + j)
            if fork_neighbor == smc_pos:
                return True
        return False
    
    def _push_off_smc_at_pos(self, pos):
        """
        Push off (dissociate) any SMC that has an anchor within smc_width of position pos.
        This handles the case where a replication fork moves through an SMC.
        Returns list of SMC indices that were pushed off.
        """
        pushed_off = []
        for i in range(self.num_smc):
            if self.loops1[i] == 0 and self.loops2[i] == 0:
                continue  # Skip inactive loops
            
            # Check if either anchor is within smc_width of the fork position
            pos1_near = self._is_pos_near_fork(self.loops1[i], pos)
            pos2_near = self._is_pos_near_fork(self.loops2[i], pos)
            
            if pos1_near or pos2_near:
                # This SMC is in the way of the fork - push it off
                self.release_around(self.loops1[i])
                self.release_around(self.loops2[i])
                self.stalled1[i] = 0
                self.stalled2[i] = 0
                self.birth(i)  # Rebirth at new location
                pushed_off.append(i)
        return pushed_off
    
    def replicate_one_step(self):
        """
        Replicate the chromosome one step.
        Pushes off SMCs that are in the way of moving forks.
        """
        # Get old fork positions before replication
        old_left_fork = self.chromo.get_left_fork_position()
        old_right_fork = self.chromo.get_right_fork_position()
        
        # Store old max_bead before replication (needed for proper re-indexing)
        old_max_bead = self.current_max_bead
        
        # Perform replication
        result = self.chromo.replicate_one_step()
        
        # Update current_max_bead (two new beads were created)
        # next_bead_id is the next available ID, so max is next_bead_id - 1
        self.current_max_bead = self.chromo.next_bead_id - 1
        
        # Get new fork positions after replication
        new_left_fork = self.chromo.get_left_fork_position()
        new_right_fork = self.chromo.get_right_fork_position()
        
        # Push off SMCs that are at the new fork positions
        pushed_off_left = self._push_off_smc_at_pos(new_left_fork)
        pushed_off_right = self._push_off_smc_at_pos(new_right_fork)
        
        # Push off SMCs that the forks passed through
        # Check positions between old and new fork positions
        if old_left_fork != new_left_fork:
            # Left fork moves backward (decreasing index, wrapping around)
            # Calculate path from old to new (backward direction)
            max_bead = self.current_max_bead
            if old_left_fork > new_left_fork:
                # Simple backward movement, no wraparound
                for pos in range(new_left_fork, old_left_fork + 1):
                    self._push_off_smc_at_pos(pos)
            else:
                # Wraparound: from old_left_fork backward to 0, then to new_left_fork
                for pos in range(old_left_fork, max_bead + 1):
                    self._push_off_smc_at_pos(pos)
                for pos in range(0, new_left_fork + 1):
                    self._push_off_smc_at_pos(pos)
        
        if old_right_fork != new_right_fork:
            # Right fork moves forward (increasing index)
            # Simple forward movement
            max_bead = self.current_max_bead
            if new_right_fork > old_right_fork:
                for pos in range(old_right_fork, new_right_fork + 1):
                    self._push_off_smc_at_pos(pos)
            else:
                # Shouldn't happen, but handle wraparound case
                for pos in range(old_right_fork, max_bead + 1):
                    self._push_off_smc_at_pos(pos)
                for pos in range(0, new_right_fork + 1):
                    self._push_off_smc_at_pos(pos)
        
        # Re-index SMCs on the right daughter to match shifted bead IDs
        # All existing right daughter beads have been shifted by +1
        # So SMC loop positions on right daughter also shift by +1
        for i in range(self.num_smc):
            if self.loops1[i] == 0 and self.loops2[i] == 0:
                continue  # Skip inactive loops
            
            # Check if anchors are on the right daughter
            if self.loops1[i] >= self.initial_length:
                # Shift by +1 to match the shifted bead ID
                self.loops1[i] = self.loops1[i] + 1
            if self.loops2[i] >= self.initial_length:
                # Shift by +1 to match the shifted bead ID
                self.loops2[i] = self.loops2[i] + 1
        
        # Rebuild occupied array from scratch (cleanest approach after re-indexing)
        # This avoids issues with wraparound behavior changing when max_bead changes
        max_beads = 2 * self.initial_length
        self.occupied = [0] * (max_beads + 1)
        
        # Re-add occupancy from all current SMC positions
        for i in range(self.num_smc):
            if self.loops1[i] != 0 or self.loops2[i] != 0:
                self.occupy_around(self.loops1[i])
                self.occupy_around(self.loops2[i])
        
        # Invalidate cache when chromosome structure changes
        self._beads_cache_invalid = True
        return result
    
    def is_fully_replicated(self):
        """Check if replication is complete."""
        return self.chromo.is_fully_replicated()

