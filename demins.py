#!/usr/bin/env python
# coding: utf-8

import pint
from io import StringIO
import datetime as dt

ureg = pint.UnitRegistry()
ureg.default_format = "0.2E"
ureg.define('ppm = 1000 gram / 1E9 gram')
ureg.define('ppb = ppm * 1E-3')
ureg.define('ppt = ppm * 1E-6')
ureg.define('gpm = gallon / minute')
ureg.define('density = gram / gallon')
ureg.define('equivalents = 1 = eq')
ureg.define('mass_conc = 1000 gram / 1000 gram')
ureg.define('i131 = mass_conc / (1e6 * 1.25E5)')
ureg.define('cs137 = mass_conc / (1e6 * 88)')
ureg.define('co58 = mass_conc / (1e6 * 3.19E4)')
ureg.define('co60 = mass_conc / (1e6 * 1.13E3)')

class Resin():
    def __init__():
        self.name = ''
        self.type = ''
        self.capacity = ''
        
class Species(object):
    def __init__(self, conc=0, name='', units='ppb', equivalent_per_mole=0, molecular_weight=0, removal_efficiency=0):
        self.name = name
        ureg_from_str = dict(ppb=ureg.ppb,
                ppm=ureg.ppm,
                i131=ureg.i131,
                cs137=ureg.cs137,
                co58=ureg.co58,
                co60=ureg.co60)
        self.equivalent_per_mole = equivalent_per_mole * ureg.eq / ureg.mole
        self.molecular_weight = molecular_weight * ureg.gram / ureg.mole
        self.removal_efficiency = removal_efficiency
        self.units = ureg_from_str[units.lower()]
        self.conc = conc * self.units
        
    def dot_e_calc(self, flow, density):
        flow = flow# * ureg.gallon / ureg.minute
        density = density * ureg.gram / ureg.gallon
        self.dot_e = self.conc.to('ppb') * self.equivalent_per_mole * flow * density / (self.molecular_weight * 1E9 * ureg.ppb)
        
    def uptake_calc(self):
        self.uptake_rate = self.dot_e * self.removal_efficiency
        

class Vessel():
    def __init__(self):
        self.anions = []
        self.cations = []
        self.cation_vol = 0
        self.anion_vol = 0
        self.start_date = ''
        self.throughput = 0
        self.service_time = 0
        self.cation_capacity_init = 0
        self.anion_capacity_init = 0
        self.cation_capacity_penalty = 0
        self.anion_capacity_penalty = 0

    def __repr__(self):
        buf = StringIO("")
        if self.anion_vol != 0:
            buf.write(f"Anion {round(self.anion_percent_exhausted().magnitude,2)}% exh with {round(self.time_to_anion_exhaust().to('days').magnitude)} days left. ")
        if self.cation_vol != 0:
            buf.write(f"Cation {round(self.cation_percent_exhausted().magnitude,2)}% exh with {round(self.time_to_cation_exhaust().to('days').magnitude)} days left.")
        return buf.getvalue()

    def initial_bed_capacity(self):
        """
        Initial resin loading capacity.
        """
        self.cation_capacity_init = (self.cation_vol.to('L') * self.cation_specific_capacity).to('eq')
        self.anion_capacity_init = (self.anion_vol.to('L') * self.anion_specific_capacity).to('eq')

    def cation_percent_exhausted(self):
        """
        Percent cation bed exhaustion.
        """
        catcap = self.cation_capacity_calc()
        try:
            self.cation_pct_exh = 100 * (1 - (catcap / self.cation_capacity_init)) #* ureg.percent
        except ZeroDivisionError:
            self.cation_pct_exh = 0 * ureg.equivalents / ureg.equivalents
        return self.cation_pct_exh

    def anion_percent_exhausted(self):
        """
        Percent anion bed exhaustion.
        """
        anicap = self.anion_capacity_calc()
        try:
            self.anion_pct_exh = 100 * (1 - (anicap / self.anion_capacity_init)) #* ureg.percent
        except ZeroDivisionError:
            self.anion_pct_exh = 0 * ureg.equivalents / ureg.equivalents
        return self.anion_pct_exh

    def time_to_cation_exhaust(self, percent_exhausted=0.7):
        """Determine time to percent bed exhaustion in seconds. """
        if self.cation_vol == 0:
            return 0 * ureg.days
        return ((percent_exhausted * self.cation_capacity_init) - (
            self.cation_capacity_init - self.cation_capacity_calc()))/self.cation_uptake_rate_calc()
    
    def time_to_anion_exhaust(self, percent_exhausted=0.7):
        """Determine time to percent bed exhaustion in seconds. """
        if self.anion_vol == 0:
            return 0 * ureg.days
        return ((percent_exhausted * self.anion_capacity_init) - (
            self.anion_capacity_init - self.anion_capacity_calc()))/self.anion_uptake_rate_calc()

    def cation_uptake_rate_calc(self):
        """ uptake rate in equivalents / minute """
        self.cation_uptake_rate = sum([sp.uptake_rate for sp in self.cations])
        return self.cation_uptake_rate

    def cation_capacity_calc(self):
        """ cation capacity in equivalents """
        t = self.service_time_calc().to('minute')
        return (self.cation_capacity_init * (1-self.cation_capacity_penalty)) - ( self.cation_uptake_rate_calc() * t)

    def anion_uptake_rate_calc(self):
        """ uptake rate in equivalents / minute """
        self.anion_uptake_rate = sum([sp.uptake_rate for sp in self.anions])
        return self.anion_uptake_rate

    def anion_capacity_calc(self):
        """ anion capacity in equivalents """
        t = self.service_time_calc().to('minute')
        return (self.anion_capacity_init * (1-self.anion_capacity_penalty)) - ( self.anion_uptake_rate_calc() * t)

    def service_time_calc(self):
        """ service time since start_date in seconds"""
        self.service_time = (dt.datetime.now() - self.start_date).total_seconds() * ureg.second
        return self.service_time

    def vessel_throughput(self):
        """ Total system volume throughput based on flow rate and time in service."""
        if self.service_time == 0:
            self.service_time_calc()
        self.throughput = self.service_time * self.flow * ureg.gal / ureg.min * ureg.min / 60 * ureg.sec
        return self.throughput

    def max_pct_exhausted(self):
        if self.cation_vol == 0:
            pct = self.anion_percent_exhausted().magnitude
        elif self.anion_vol == 0:
            pct = self.cation_percent_exhausted().magnitude
        else:
            pct = max([self.cation_percent_exhausted().magnitude,self.anion_percent_exhausted().magnitude])
        self.pct = round(pct)
        return self.pct

