#!/usr/bin/env python3

import numpy as np
import sys
import os


btu_per_hour_to_tons = 1.0/11917
watts_to_btu_per_hour = 3.41
joules_to_btu = 0.000947817


if (len(sys.argv) == 2):
  fname = sys.argv[1]
else:
  print("Usage: " + sys.argv[0] + " fname")
  exit(1)

def getColumnNames(fname):
  file = open(fname, 'r')
  line = file.readline()

  names = line.split(' ')
  names[-1] = names[-1].rstrip('\n')
  return names


def getColumnIndicesFromName(column_names, column_names_to_get):
  indices = np.zeros(0, dtype=int)
  for i in range(len(column_names_to_get)):
    found_flag = False
    for j in range(len(column_names)):
      if column_names_to_get[i] == column_names[j]:
        indices = np.append(indices, j)
        found_flag = True
        break

    if (not found_flag):
      raise ValueError("unable to find column named ", column_names_to_get[i])

  return indices


def computeTotalHVACFluxes(col_names, data):

  indices = getColumnIndicesFromName(col_names, ["time", "hvac_flux", "exterior_air_temp"])
  time = data[:, indices[0]]
  hvac_flux = data[:, indices[1]]
  temp  = data[:, indices[2]]
  heating_flux = np.maximum(hvac_flux, 0)
  cooling_flux = np.minimum(hvac_flux, 0)

  net_heating_flux = np.trapz(heating_flux, time)
  net_cooling_flux = np.trapz(cooling_flux, time)

  net_heating_flux_btu = net_heating_flux * joules_to_btu/(1000*1000)
  net_cooling_flux_btu = net_cooling_flux * joules_to_btu/(1000*1000)

  max_heating_flux = np.max(heating_flux)
  max_cooling_flux = np.min(cooling_flux)

  max_t = np.max(temp)
  min_t = np.min(temp)

  max_t_f = (max_t - 273.15) * 9.0/5 + 32
  min_t_f = (min_t - 273.15) * 9.0/5 + 32

  max_heating_flux_tons = max_heating_flux * watts_to_btu_per_hour * btu_per_hour_to_tons
  max_cooling_flux_tons = max_cooling_flux * watts_to_btu_per_hour * btu_per_hour_to_tons
  print("  net heating flux = ", net_heating_flux/(1000*3600), " kWh, or ", net_heating_flux_btu, " mBTU")
  print("  max heating flux = ", max_heating_flux, " W, or ", max_heating_flux_tons, " tons")
  print("  net cooling flux = ", net_cooling_flux/(1000*3600), " kWh, or ", net_cooling_flux_btu, " mBTU")
  print("  max cooling flux = ", max_cooling_flux, " W, or ", max_cooling_flux_tons, " tons")
  print("  temperature range = ", min_t_f, ", ", max_t_f, " F")


column_names = getColumnNames(fname)
data = np.loadtxt(fname, skiprows=1)

computeTotalHVACFluxes(column_names, data)

