"""
Formate Beads Release Calculator - Core Logic

This module contains the core calculation logic for formate bead experiments.
Each bead releases formate over multiple days with different release rates.
Fast beads: 80 mg total formate
Slow beads: 60 mg total formate
"""

# =============================================================================
# CONFIGURATION - EDIT THE RELEASE RATES HERE
# =============================================================================

# Release rate for M07 beads (mmol/day for each day)
# These values represent the daily release rate from ONE bead
M07_BEAD_RELEASE = {
    1: 0.300,  # Day 1: 0.300 mmol/day
    2: 0.250,  # Day 2: 0.250 mmol/day
    3: 0.200,  # Day 3: 0.200 mmol/day
    4: 0.150,  # Day 4: 0.150 mmol/day
    5: 0.120,  # Day 5: 0.120 mmol/day
    6: 0.090,  # Day 6: 0.090 mmol/day
    7: 0.066,  # Day 7: 0.066 mmol/day
}

# Release rate for M03 beads (mmol/day for each day)
# These values represent the daily release rate from ONE bead
M03_BEAD_RELEASE = {
    1: 0.150,  # Day 1: 0.150 mmol/day
    2: 0.140,  # Day 2: 0.140 mmol/day
    3: 0.130,  # Day 3: 0.130 mmol/day
    4: 0.120,  # Day 4: 0.120 mmol/day
    5: 0.110,  # Day 5: 0.110 mmol/day
    6: 0.100,  # Day 6: 0.100 mmol/day
    7: 0.080,  # Day 7: 0.080 mmol/day
}

# =============================================================================
# CONSTANTS
# =============================================================================

FORMATE_MW = 68  # mg/mmol
TOTAL_FORMATE_M07_BEAD = 80  # mg per M07 bead
TOTAL_FORMATE_M03_BEAD = 60  # mg per M03 bead
TOTAL_MMOL_M07_BEAD = TOTAL_FORMATE_M07_BEAD / FORMATE_MW  # 1.176 mmol
TOTAL_MMOL_M03_BEAD = TOTAL_FORMATE_M03_BEAD / FORMATE_MW  # 0.882 mmol


# =============================================================================
# BEAD CLASS
# =============================================================================

class Bead:
    """
    Represents a single formate-releasing bead.
    
    Tracks when it was added, how much it has released, and when it's depleted.
    """
    
    def __init__(self, day_added, bead_type, release_profile):
        """
        Initialize a bead.
        
        Parameters:
        -----------
        day_added : int
            Day when this bead was added to the experiment
        bead_type : str
            'M07' or 'M03'
        release_profile : dict
            Dictionary mapping day number to release rate (mmol/day)
        """
        self.day_added = day_added
        self.bead_type = bead_type
        self.release_profile = release_profile
        # Set capacity based on bead type
        if bead_type == 'M07':
            self.total_capacity = TOTAL_MMOL_M07_BEAD  # 1.176 mmol
        else:
            self.total_capacity = TOTAL_MMOL_M03_BEAD  # 0.882 mmol
        self.total_released = 0.0  # Track cumulative release
        self.depleted = False
    
    def get_release_rate(self, current_day):
        """
        Get the release rate for a specific day.
        
        Parameters:
        -----------
        current_day : int
            Current day of experiment
        
        Returns:
        --------
        float : Release rate in mmol/day (0 if bead is depleted or outside release window)
        """
        if self.depleted:
            return 0.0
        
        days_since_added = current_day - self.day_added + 1
        
        # Check if within release profile range
        if days_since_added < 1 or days_since_added > len(self.release_profile):
            return 0.0
        
        release_rate = self.release_profile.get(days_since_added, 0.0)
        
        # Check if this would exceed capacity
        if self.total_released + release_rate > self.total_capacity:
            # Release only remaining capacity
            remaining = self.total_capacity - self.total_released
            self.total_released = self.total_capacity
            self.depleted = True
            return remaining
        
        self.total_released += release_rate
        return release_rate
    
    def is_depleted(self):
        """Check if bead has released all its formate."""
        return self.depleted
    
    def __repr__(self):
        return (f"Bead(type={self.bead_type}, added_day={self.day_added}, "
                f"released={self.total_released:.3f}/{self.total_capacity:.3f} mmol)")


# =============================================================================
# EXPERIMENT MANAGER CLASS
# =============================================================================

class ExperimentManager:
    """
    Manages the formate release experiment across multiple days.
    Supports mixing M07 and M03 beads on the same day.
    Tracks cumulative formate concentration in the medium.
    """
    
    def __init__(self, volume_ml):
        """
        Initialize experiment manager.
        
        Parameters:
        -----------
        volume_ml : float
            Volume of the experimental medium in mL
        """
        self.volume_ml = volume_ml
        self.volume_L = volume_ml / 1000.0  # Convert to liters for mM calculation
        
        self.active_beads = []  # List of Bead objects
        self.schedule = {}  # Day -> {'M07': count, 'M03': count}
        self.cumulative_formate_mmol = 0.0  # Track cumulative formate in medium
    
    def calculate_beads_needed(self, desired_cumulative_concentration_mM):
        """
        Calculate how many beads to add each day to reach desired CUMULATIVE concentration.
        Uses a combination of M07 and M03 beads for optimal flexibility.
        
        Parameters:
        -----------
        desired_cumulative_concentration_mM : dict
            Dictionary with day (1-7) as key and desired CUMULATIVE concentration (mM) as value
            Example: {1: 5.0, 2: 10.0, 3: 15.0, ...} - concentration builds up over time
        
        Returns:
        --------
        dict : Schedule of beads to add each day
        """
        print(f"\n{'='*80}")
        print(f"FORMATE BEAD CALCULATION - MIXED M07 & M03 BEADS")
        print(f"CUMULATIVE CONCENTRATION MODE")
        print(f"{'='*80}\n")
        print(f"Experimental Volume: {self.volume_ml} mL ({self.volume_L} L)")
        print(f"\nDesired Cumulative Concentrations:")
        for day in sorted(desired_cumulative_concentration_mM.keys()):
            conc_mM = desired_cumulative_concentration_mM[day]
            mmol_needed = conc_mM * self.volume_L
            print(f"  Day {day}: {conc_mM:.2f} mM cumulative ({mmol_needed:.3f} mmol total)")
        print(f"\n{'='*80}\n")
        
        for day in range(1, 8):
            if day not in desired_cumulative_concentration_mM:
                print(f"Day {day}: No target specified - skipping\n")
                continue
            
            desired_cumulative_conc = desired_cumulative_concentration_mM[day]
            desired_cumulative_mmol = desired_cumulative_conc * self.volume_L  # Convert mM to mmol total
            
            # Calculate today's release from all active beads
            todays_release_mmol = 0
            active_count = 0
            depleted_count = 0
            
            for bead in self.active_beads:
                release = bead.get_release_rate(day)
                todays_release_mmol += release
                if not bead.is_depleted():
                    active_count += 1
                else:
                    depleted_count += 1
            
            # Update cumulative formate
            self.cumulative_formate_mmol += todays_release_mmol
            current_cumulative_conc = self.cumulative_formate_mmol / self.volume_L if self.volume_L > 0 else 0
            
            # Calculate additional release needed to reach target cumulative
            additional_mmol_needed = desired_cumulative_mmol - self.cumulative_formate_mmol
            
            # Initialize bead counts
            m07_beads_needed = 0
            m03_beads_needed = 0
            
            print(f"Day {day}:")
            print(f"  Target cumulative:      {desired_cumulative_conc:.2f} mM ({desired_cumulative_mmol:.3f} mmol)")
            print(f"  Today's release:        {todays_release_mmol / self.volume_L:.2f} mM ({todays_release_mmol:.3f} mmol)")
            print(f"  Current cumulative:     {current_cumulative_conc:.2f} mM ({self.cumulative_formate_mmol:.3f} mmol)")
            print(f"  Active beads:           {active_count} beads")
            if depleted_count > 0:
                print(f"  Depleted beads:         {depleted_count} beads")
            
            if additional_mmol_needed <= 0:
                print(f"  → No beads needed (already at or above target)")
                # Calculate actual vs target
                deviation_mmol = self.cumulative_formate_mmol - desired_cumulative_mmol
                deviation_conc = current_cumulative_conc - desired_cumulative_conc
                if abs(deviation_mmol) > 0.001:
                    print(f"  Deviation: {deviation_conc:+.3f} mM ({deviation_mmol:+.3f} mmol)")
            else:
                additional_conc = additional_mmol_needed / self.volume_L
                print(f"  Additional needed:      {additional_conc:.2f} mM ({additional_mmol_needed:.3f} mmol)")
                
                # Use combination: prioritize M03 beads for sustained release
                # Strategy: use M03 beads for baseline, M07 beads for extra boost
                m07_day1_release = M07_BEAD_RELEASE[1]
                m03_day1_release = M03_BEAD_RELEASE[1]
                
                # Calculate ideal combination
                m03_beads_float = additional_mmol_needed / m03_day1_release
                
                # Round to whole beads and optimize combination
                if m03_beads_float > 3:
                    # Use 2 M03 beads and make up difference with M07
                    m03_beads_needed = 2
                    remaining = additional_mmol_needed - (m03_beads_needed * m03_day1_release)
                    m07_beads_needed = int(round(remaining / m07_day1_release))
                else:
                    # Round M03 beads
                    m03_beads_needed = int(round(m03_beads_float))
                    m07_beads_needed = 0
                
                # Calculate actual cumulative with whole beads
                actual_additional = (m07_beads_needed * m07_day1_release + 
                                   m03_beads_needed * m03_day1_release)
                actual_cumulative = self.cumulative_formate_mmol + actual_additional
                actual_cumulative_conc = actual_cumulative / self.volume_L
                
                # Update cumulative formate with new beads
                self.cumulative_formate_mmol = actual_cumulative
                
                # Calculate deviation
                deviation_mmol = actual_cumulative - desired_cumulative_mmol
                deviation_conc = actual_cumulative_conc - desired_cumulative_conc
                
                print(f"  → ADD {m07_beads_needed} M07 beads + {m03_beads_needed} M03 beads")
                print(f"  Actual cumulative: {actual_cumulative_conc:.3f} mM ({actual_cumulative:.3f} mmol)")
                print(f"  Deviation: {deviation_conc:+.3f} mM ({deviation_mmol:+.3f} mmol)")
                
                # Add new beads
                for _ in range(m07_beads_needed):
                    self.active_beads.append(Bead(day, 'M07', M07_BEAD_RELEASE))
                for _ in range(m03_beads_needed):
                    self.active_beads.append(Bead(day, 'M03', M03_BEAD_RELEASE))
            
            self.schedule[day] = {'M07': m07_beads_needed, 'M03': m03_beads_needed}
            print()
        
        self._print_summary()
        return self.schedule
    
    def _print_summary(self):
        """Print experiment summary."""
        print(f"{'='*80}")
        print(f"SUMMARY")
        print(f"{'='*80}\n")
        
        total_m07 = sum(day['M07'] for day in self.schedule.values())
        total_m03 = sum(day['M03'] for day in self.schedule.values())
        total_beads = total_m07 + total_m03
        
        total_formate_mmol = self.cumulative_formate_mmol
        total_formate_mg = total_formate_mmol * FORMATE_MW
        
        print(f"Total M07 beads added:    {total_m07} beads")
        print(f"Total M03 beads added:    {total_m03} beads")
        print(f"Total beads:              {total_beads} beads")
        print(f"Final cumulative formate: {total_formate_mmol:.3f} mmol ({total_formate_mg:.2f} mg)")
        print(f"Final concentration:      {(total_formate_mmol / self.volume_L):.2f} mM")
        print(f"\nBead efficiency check:")
        
        m07_beads = [b for b in self.active_beads if b.bead_type == 'M07']
        m03_beads = [b for b in self.active_beads if b.bead_type == 'M03']
        
        print(f"  M07 beads ({len(m07_beads)} total):")
        if m07_beads:
            fully_used = sum(1 for b in m07_beads if b.total_released > 0.9 * b.total_capacity)
            partially_used = sum(1 for b in m07_beads if 0.1 * b.total_capacity < b.total_released <= 0.9 * b.total_capacity)
            barely_used = sum(1 for b in m07_beads if b.total_released <= 0.1 * b.total_capacity)
            print(f"    Fully utilized (>90%):    {fully_used} beads")
            print(f"    Partially used (10-90%):  {partially_used} beads")
            print(f"    Barely used (<10%):       {barely_used} beads")
        
        print(f"  M03 beads ({len(m03_beads)} total):")
        if m03_beads:
            fully_used = sum(1 for b in m03_beads if b.total_released > 0.9 * b.total_capacity)
            partially_used = sum(1 for b in m03_beads if 0.1 * b.total_capacity < b.total_released <= 0.9 * b.total_capacity)
            barely_used = sum(1 for b in m03_beads if b.total_released <= 0.1 * b.total_capacity)
            print(f"    Fully utilized (>90%):    {fully_used} beads")
            print(f"    Partially used (10-90%):  {partially_used} beads")
            print(f"    Barely used (<10%):       {barely_used} beads")
        print()
