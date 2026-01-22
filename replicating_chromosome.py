"""
ReplicatingChromosome: A circular chromosome that can replicate bidirectionally from an origin.

Initially, the chromosome has N beads (0 to N-1) in a circular arrangement.
Replication starts at the origin (bead 0) and proceeds bidirectionally.
"""

import os


class ReplicatingChromosome:
    """
    Represents a circular chromosome that can replicate bidirectionally from an origin.
    
    Initially, the chromosome has N beads (0 to N-1) in a circular arrangement.
    Replication starts at the origin (bead 0) and proceeds bidirectionally.
    """
    
    def __init__(self, initial_length, coords=None, box_size=(200.0, 200.0, 200.0)):
        """
        Initialize an unreplicated circular chromosome.
        
        Args:
            initial_length: Number of beads in the initial chromosome (e.g., 5000)
            coords: Optional dict mapping bead_id -> (x, y, z) tuple. If None, coordinates not initialized.
            box_size: Tuple of (x, y, z) box dimensions (default: (200.0, 200.0, 200.0))
        """
        self.initial_length = initial_length
        self.next_bead_id = initial_length  # Next available bead ID for new beads
        self.replicated_left = 0   # Number of beads replicated to the left (backwards from origin)
        self.replicated_right = 0  # Number of beads replicated to the right (forwards from origin)
        
        # Store bonds as a set of tuples (bead1, bead2) where bead1 < bead2
        self.bonds = set()
        
        # Initialize initial circular chromosome bonds
        for i in range(initial_length):
            next_bead = (i + 1) % initial_length
            if i < next_bead:
                self.bonds.add((i, next_bead))
            else:
                self.bonds.add((next_bead, i))
        
        # Coordinate storage: dict mapping bead_id -> (x, y, z) tuple
        self.coords = coords if coords is not None else {}
        self.box_size = box_size
        
        # Train track displacement: how much to displace daughter beads from mother
        self.daughter_displacement = 0.1  # Can be adjusted
    
    def get_left_fork_position(self):
        """Get the bead index of the left replication fork."""
        # After replicating n beads to the left, the fork is at position (initial_length - n)
        # For example, after 1 step, fork is at 4998 (which is 5000 - 1 - 1)
        return (self.initial_length - self.replicated_left - 1) % self.initial_length
    
    def get_right_fork_position(self):
        """Get the bead index of the right replication fork."""
        # After replicating n beads to the right, the fork is at position n
        # For example, after 1 step, fork is at 1
        return self.replicated_right
    
    def replicate_one_step(self):
        """
        Replicate one bead in both directions from the replication forks.
        
        Returns:
            tuple: (left_new_bead, right_new_bead) - IDs of the two new beads created
        """
        if self.is_fully_replicated():
            raise ValueError("Chromosome is already fully replicated")
        
        if self.replicated_left == 0 and self.replicated_right == 0:
            current_left_fork = 0
            current_right_fork = 0
        else:
            current_left_fork = self.get_left_fork_position()
            current_right_fork = self.get_right_fork_position()
        
        if self.replicated_left == 0 and self.replicated_right == 0:
            new_right_fork = 1
            new_left_fork = (self.initial_length - 2) % self.initial_length
        else:
            new_right_fork = (current_right_fork + 1) % self.initial_length
            new_left_fork = (current_left_fork - 1) % self.initial_length
        
        existing_right_beads = [b for b in self.get_all_beads() if b >= self.initial_length]
        existing_right_beads.sort()
        
        left_new_bead = self.initial_length
        
        if existing_right_beads:
            right_new_bead = max(existing_right_beads) + 2
        else:
            right_new_bead = self.initial_length + 1
        
        bead_id_mapping = {}
        temp_coords_storage = {}
        
        for old_id in existing_right_beads:
            new_id = old_id + 1
            bead_id_mapping[old_id] = new_id
            if old_id in self.coords:
                temp_coords_storage[new_id] = self.coords[old_id]
        
        for old_id in reversed(existing_right_beads):
            new_id = bead_id_mapping[old_id]
            if new_id in temp_coords_storage:
                self.coords[new_id] = temp_coords_storage[new_id]
        
        for old_id in existing_right_beads:
            new_id = bead_id_mapping[old_id]
            if (old_id != new_id and 
                new_id in temp_coords_storage and 
                new_id in self.coords):
                if old_id in self.coords:
                    is_also_destination = old_id in bead_id_mapping.values()
                    if not is_also_destination:
                        del self.coords[old_id]
        
        self.next_bead_id = right_new_bead + 1
        
        is_first_step = (self.replicated_left == 0 and self.replicated_right == 0)
        will_be_fully_replicated = (self.replicated_left + self.replicated_right + 2) >= self.initial_length
        
        if is_first_step:
            self.bonds.add(tuple(sorted([left_new_bead, right_new_bead])))
            self.bonds.add(tuple(sorted([new_left_fork, left_new_bead])))
            self.bonds.add(tuple(sorted([right_new_bead, new_right_fork])))
        elif not will_be_fully_replicated:
            current_rightmost_right_daughter = right_new_bead - 2
            first_right_daughter_bead = self.initial_length
            
            current_right_fork_bond = tuple(sorted([current_rightmost_right_daughter, current_right_fork]))
            self.bonds.discard(current_right_fork_bond)
            current_left_fork_bond = tuple(sorted([current_left_fork, first_right_daughter_bead]))
            self.bonds.discard(current_left_fork_bond)
            
            self.bonds.add(tuple(sorted([right_new_bead-1, right_new_bead])))
            self.bonds.add(tuple(sorted([current_rightmost_right_daughter, right_new_bead-1])))
            self.bonds.add(tuple(sorted([new_left_fork, first_right_daughter_bead])))
            self.bonds.add(tuple(sorted([right_new_bead, new_right_fork])))
        else:
            print("Final replication step")
            current_rightmost_right_daughter = right_new_bead - 2
            first_right_daughter_bead = self.initial_length
            
            current_right_fork_bond = tuple(sorted([current_rightmost_right_daughter, current_right_fork]))
            self.bonds.discard(current_right_fork_bond)         
            current_left_fork_bond = tuple(sorted([current_left_fork, first_right_daughter_bead]))
            self.bonds.discard(current_left_fork_bond)
            
            self.bonds.add(tuple(sorted([right_new_bead-1, right_new_bead])))
            self.bonds.add(tuple(sorted([current_rightmost_right_daughter, right_new_bead-1])))
            self.bonds.add(tuple(sorted([right_new_bead, first_right_daughter_bead])))
        
        self.replicated_left += 1
        self.replicated_right += 1
        
        if new_left_fork in self.coords:
            x, y, z = self.coords[new_left_fork]
            self.coords[left_new_bead] = (x, y, z + self.daughter_displacement)
        else:
            self.coords[left_new_bead] = (0.0, 0.0, self.daughter_displacement)
        
        if new_right_fork in self.coords:
            x, y, z = self.coords[new_right_fork]
            self.coords[right_new_bead] = (x, y, z + self.daughter_displacement)
        else:
            self.coords[right_new_bead] = (0.0, 0.0, self.daughter_displacement)
        
        return (left_new_bead, right_new_bead)
    
    def is_fully_replicated(self):
        """
        Check if replication is complete (two separate circular chromosomes).
        
        Returns:
            bool: True if fully replicated
        """
        # Replication is complete when we've replicated half the chromosome in each direction
        return (self.replicated_left + self.replicated_right) >= self.initial_length
    
    def get_all_beads(self):
        """Get a list of all bead IDs."""
        # Original beads: 0 to initial_length-1
        # New beads: initial_length to next_bead_id-1
        total_beads = self.next_bead_id
        return list(range(total_beads))
    
    def get_bonds(self):
        """Get all bonds as a list of tuples."""
        return sorted(list(self.bonds))
    
    def load_coords_from_file(self, filename, vertex_offset=0):
        """
        Load coordinates from a dna.q format file.
        
        Args:
            filename: Path to dna.q file
            vertex_offset: Offset to subtract from vertex IDs to get internal bead IDs (default: 0)
                          This accounts for membrane vertices that come before DNA beads
        """
        with open(filename, 'r') as f:
            lines = f.readlines()
        
        # Parse box size (first line)
        if len(lines) > 0:
            box_values = list(map(float, lines[0].split()))
            if len(box_values) >= 3:
                self.box_size = tuple(box_values[:3])
        
        # Parse number of beads (second line)
        num_beads = int(lines[1].strip()) if len(lines) > 1 else 0
        
        # Parse coordinates (lines 2 to num_beads+1)
        self.coords = {}
        for i in range(num_beads):
            if i + 2 < len(lines):
                parts = lines[i + 2].split()
                if len(parts) >= 4:
                    vertex_id = int(parts[0])
                    x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                    # Convert vertex ID to internal bead ID by subtracting offset
                    bead_id = vertex_id - vertex_offset
                    # Only store coordinates for DNA beads (bead_id >= 0)
                    if bead_id >= 0:
                        self.coords[bead_id] = (x, y, z)
        
        # Ensure we have coordinates for all initial beads
        for i in range(self.initial_length):
            if i not in self.coords:
                # Default to origin if not specified
                self.coords[i] = (0.0, 0.0, 0.0)
    
    def write_coords_to_file(self, filename, vertex_offset=0):
        """
        Write coordinates to a dna.q format file.
        
        Args:
            filename: Path to output dna.q file
            vertex_offset: Offset to add to all DNA bead IDs (default: 0)
                          This accounts for membrane vertices that come before DNA beads
        """
        all_beads = sorted(self.get_all_beads())
        
        with open(filename, 'w') as f:
            # Write box size
            f.write(f"{self.box_size[0]:.10f}   {self.box_size[1]:.10f}   {self.box_size[2]:.10f}\n")
            
            # Write number of beads
            f.write(f"{len(all_beads)}\n")
            
            # Write coordinates for each bead (with offset)
            for bead_id in all_beads:
                vertex_id = bead_id + vertex_offset
                if bead_id in self.coords:
                    x, y, z = self.coords[bead_id]
                    f.write(f"{vertex_id}   {x:.10f}   {y:.10f}   {z:.10f}  0\n")
                else:
                    # Default to origin if coordinate not set
                    f.write(f"{vertex_id}   0.0000000000   0.0000000000   0.0000000000  0\n")
            
            # Write terminator
            f.write("0\n")
            # Ensure file is flushed to disk
            f.flush()
            os.fsync(f.fileno())
    
    def visualize_state(self):
        """Print current state for debugging."""
        print(f"Initial length: {self.initial_length}")
        print(f"Replicated left: {self.replicated_left}, Replicated right: {self.replicated_right}")
        print(f"Left fork at: {self.get_left_fork_position()}, Right fork at: {self.get_right_fork_position()}")
        print(f"Total beads: {len(self.get_all_beads())}")
        print(f"Total bonds: {len(self.bonds)}")
        print(f"Bonds: {self.get_bonds()[:10]}..." if len(self.bonds) > 10 else f"Bonds: {self.get_bonds()}")

