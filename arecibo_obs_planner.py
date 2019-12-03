#!/usr/bin/env python3
from astropy.coordinates import EarthLocation as E
from astropy import units as u
from astropy.table import Table
from astropy.time import Time
from astropy.coordinates import SkyCoord, Angle, AltAz

from astroplan import AltitudeConstraint
from astroplan import Observer, FixedTarget
from astroplan import observability_table
from astroplan import time_grid_from_range
from astropy.time import Time, TimezoneInfo

import pandas as pd
import argparse
import numpy as np
from tabulate import tabulate

from datetime import datetime

def make_obs_table(obs_info, source_list, save=True, min_uptime=10):

    tab = Table.read(source_list, format="ascii.csv")
    targets = [FixedTarget(coord=SkyCoord(ra=ra, dec=dec, unit=(u.hourangle, u.deg)), name=source)
               for source, ra, dec, _ in tab]

    
    # Some arecibo specific settings
    arecibo_site = E.from_geocentric(x=2390490.0, y = -5564764.0, z=1994727.0, unit=u.meter)
    arecibo = Observer(location=arecibo_site, name="AO")
    constraints = [AltitudeConstraint((90-19.7)*u.deg, (90-1.06)*u.deg)]
    min_alt = (90-19.7)*u.deg
    max_alt = (90-1.06)*u.deg
    utc_offset = 4*u.hour
    ast_to_utc_offset = TimezoneInfo(utc_offset=utc_offset)

    months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

    tstarts = []
    tends = []
    for index, row in obs_info.iterrows():
        obs_month = row['Month']
        for i, month in enumerate(months):
            if month == obs_month:
                break
        i = i+1
        if row['End_time(AST)'] == '24:00':
            row['End_time(AST)'] = '23:59'
        tstarts.append(f"{row['Year']}-{i:02d}-{row['Start_date']} {row['Start_time(AST)']}")
        tends.append(f"{row['Year']}-{i:02d}-{row['End_date']} {row['End_time(AST)']}")
    
    ots = []
    for st, et in zip(tstarts, tends):     

        #setting up for altitude constraints
        tstart_obs = Time(st) + utc_offset
        tobs_len = (Time(et) - Time(st)).sec/(60*60) #hours
        delta_t = np.linspace(0, tobs_len, 100)*u.hour
        frame_obs = AltAz(obstime=tstart_obs+delta_t,
                                  location=arecibo_site)

        #setting up for astroplan observability table
        st_utc = str(Time(st).to_datetime(timezone=ast_to_utc_offset))[:-6]
        et_utc = str(Time(et).to_datetime(timezone=ast_to_utc_offset))[:-6]

        time_range = Time([st_utc, et_utc], scale='utc')
                     
        print(f'Making observability table for {st}')
        ot = observability_table(constraints, arecibo, targets, times=time_grid_from_range(time_range, 5*u.min)).to_pandas()

        ot['ever obs'] = ot['ever observable']

        ot = ot.drop(['fraction of time observable', 'ever observable', 'always observable'], axis=1)
        ot_obs = ot[ot['ever obs']]

        ra = []
        dec = []
        set_time = []
        up_time = []
        rise_time = []
        tstart_observable = []
        tend_observable = []
                     
        print(f'Processing observability table for {st}')
        for i,row in ot_obs.iterrows():
            ra.append(targets[i].ra.deg)
            dec.append(targets[i].dec.deg)

            # manual altitude constraints 
            target = SkyCoord(ra=targets[i].ra, dec=targets[i].dec)
            frame_obs_altaz = target.transform_to(frame_obs)

            times_observable = tstart_obs + delta_t[(frame_obs_altaz.alt > min_alt) & (frame_obs_altaz.alt < max_alt)]

            tstart_observable.append(times_observable.min() - utc_offset)

            if Time(times_observable.max()) < Time(et_utc):
                tend_observable.append(times_observable.max() - utc_offset)
                up_time.append((Time(times_observable.max()) - Time(times_observable.min())).sec/60)
            else:
                tend_observable.append(Time(et_utc) - utc_offset)
                up_time.append((Time(et_utc) - Time(times_observable.min())).sec/60)

        ot_obs['RA'] = ra
        ot_obs['DEC'] = dec
        ot_obs['uptime (min)'] = up_time

        ot_obs['tstart_obs (AST)'] = tstart_observable 
        ot_obs['tend_obs (AST)'] = tend_observable

        ot_obs = ot_obs[ot_obs['uptime (min)'] > min_uptime]
        ot_obs = ot_obs.drop(['ever obs'], axis=1)
        ots.append(ot_obs)

        if save:
            fname = str(Time(st_utc) - utc_offset)
            f = open(fname.replace(' ', '_')+'.tab', 'w')
            f.write(f'Observations from {str(Time(st_utc) - utc_offset)} to {str(Time(et_utc) - utc_offset)}\n')
            f.write(f'Length of observations: {tobs_len:.2f}hr \n')
            f.write(tabulate(ot_obs, headers='keys', tablefmt='psql', showindex="never"))
            f.close()
        
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate observation tables for arecibo observations",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-l', '--obs_link', type=str, help='NAIC link of the proposal', required=True)
    parser.add_argument('-s', '--source_list', type=str, help='Location of the source list', required=True)
    parser.add_argument('-mu', '--min_uptime', type=int, help='Minimum Uptime required (min)', default=10)
                     
    values = parser.parse_args()
    df = pd.read_html(values.obs_link)[0]
    make_obs_table(obs_info = df, source_list = values.source_list, save=True, min_uptime=values.min_uptime)
