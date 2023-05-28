#!/usr/bin/env python3

import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('AGG')
import numpy as np
import sys


btu_per_house_to_tons = 1.0/11917
watts_to_btu_per_hour = 0.2931
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

def removeTransient(data):
  return data #data[100:-1, :]


def copyColumn(data, col_names, col_name, factor, new_name):
  indices = getColumnIndicesFromName(col_names, [col_name]);
  new_data = np.zeros((data.shape[0], 1));
  new_data[:, 0] = factor * data[:, indices[0]];
  return np.append(data, new_data, axis=1), np.append(col_names, new_name);


def sumColumns(data, col_names, col_names_to_sum, new_name):
  indices = getColumnIndicesFromName(col_names, col_names_to_sum)
  new_data = np.zeros((data.shape[0], 1))
  new_data[:, 0] = np.sum(data[:, indices], axis=1)
  return np.append(data, new_data, axis=1), np.append(col_names, new_name);


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


def convertTemperatures(data, col_names):
  temp_col_names = ["interior_air_temp", "exterior_air_temp"]
  new_col_names = ["interior_air_temp_F", "exterior_air_temp_F"]
  col_idxs = getColumnIndicesFromName(col_names, temp_col_names)

  for i in range(len(temp_col_names)):
    old_data = data[:, col_idxs[i]];
    new_data = np.zeros((data.shape[0], 1))
    new_data[:, 0] = ((old_data - 273.15) * 9.0/5) + 32

    data = np.append(data, new_data, axis=1)
    col_names = np.append(col_names, new_col_names[i])

  return data, col_names


def convertFluxesToTons(data, col_names):

  temp_col_names = ["hvac_flux", "interior_wall_flux", "exterior_wall_flux"]
  new_col_names = ["hvac_flux_tons", "interior_wall_flux_tons", "exterior_wall_flux_tons"]
  col_idxs = getColumnIndicesFromName(col_names, temp_col_names)

  for i in range(len(temp_col_names)):
    old_data = data[:, col_idxs[i]];
    new_data = np.zeros((data.shape[0], 1))
    new_data[:, 0] = old_data * watts_to_btu_per_hour * btu_per_house_to_tons

    data = np.append(data, new_data, axis=1)
    col_names = np.append(col_names, new_col_names[i])

  return data, col_names


def plotColumns(column_names, column_names_to_plot, data, xlabel, ylabel, fig, ax):
  return plotColumnsWithLabels(column_names, column_names_to_plot, data, xlabel, ylabel, fig, ax, column_names_to_plot)

def plotColumnsWithLabels(column_names, column_names_to_plot, data, xlabel, ylabel, fig, ax, labels):

  indices = getColumnIndicesFromName(column_names, column_names_to_plot)
  for i in range(1, len(indices)):
    print("plotting ", column_names_to_plot[i])
    ax.plot(data[:, indices[0]], data[:, indices[i]], label=labels[i])

  ax.set_xlabel(xlabel)
  ax.set_ylabel(ylabel)
  ymin, ymax = ax.get_ylim()
  #ax.set_yticks(np.linspace(ymin, ymax, 5))
  ax.set_yticks(getAxisLabels(ymin, ymax, 8))
  #ax.minorticks_on()
  ax.grid(True, which='both')

def getAxisLabels(ymin, ymax, npts):

  if ymax < 0 or ymin > 0:
    return np.linspace(ymin, ymax, npts);
  else:
    pts = np.linspace(ymin, ymax, npts)
    delta = ymax - ymin;
    delta_plus = np.abs(ymax)
    delta_minus = np.abs(ymin)/delta

    if (delta_plus > 0.1 and delta_minus > 0.1):
      minidx = -1
      mindist = sys.float_info.max
      for i in range(1, npts-1):
        if np.abs(pts[i]) < mindist:
          minidx = i;
          mindist = np.abs(pts[i])

      pts[minidx] = 0

    return pts;
      

def makeLegend(fig, axes):

  for ax in axes:
    # make space for the legend so it isn't on top
    # of the data
    xmin, xmax = ax.get_xlim()
    new_xmax = xmax + (xmax - xmin)*0.50
    ax.set_xlim(xmin, new_xmax)

    ax.legend(loc=1)

def computeTotalHVACFluxes(col_names, data):

  indices = getColumnIndicesFromName(col_names, ["time", "hvac_flux"])
  time = data[:, indices[0]]
  hvac_flux = data[:, indices[1]]
  heating_flux = np.maximum(hvac_flux, 0)
  cooling_flux = np.minimum(hvac_flux, 0)

  net_heating_flux = np.trapz(heating_flux, time)
  net_cooling_flux = np.trapz(cooling_flux, time)

  net_heating_flux_btu = net_heating_flux * joules_to_btu/(1000*1000)
  net_cooling_flux_btu = net_cooling_flux * joules_to_btu/(1000*1000)
  print("net heating flux = ", net_heating_flux/(1000*3600), " kWh, or ", net_heating_flux_btu, " mBTU")
  print("net cooling flux = ", net_cooling_flux/(1000*3600), " kWh, or ", net_cooling_flux_btu, " mBTU")

def writeResults(column_names, data, filename):
  header=""
  for name in column_names:
    header = header + name + " "

  print("saving to ", filename)
  np.savetxt(filename, data, header=header, comments='');
    



column_names = getColumnNames(fname)
print("names = ", column_names)
print ("about to load data")
data = np.loadtxt(fname, skiprows=1)
data = removeTransient(data)

interior_flux_names = ["east_interior_wall_flux", "north_interior_wall_flux", "west_interior_wall_flux", "south_interior_wall_flux"];
data, column_names = sumColumns(data, column_names, interior_flux_names, "interior_wall_flux");

exterior_flux_names = ["east_exterior_wall_flux", "north_exterior_wall_flux", "west_exterior_wall_flux", "south_exterior_wall_flux"];
data, column_names = sumColumns(data, column_names, exterior_flux_names, "exterior_wall_flux");

exterior_flux_names = ["east_exterior_wall_flux_cond", "north_exterior_wall_flux_cond", "west_exterior_wall_flux_cond", "south_exterior_wall_flux_cond"];
data, column_names = sumColumns(data, column_names, exterior_flux_names, "exterior_wall_flux_cond");

exterior_flux_names = ["east_exterior_wall_flux_solar_rad", "north_exterior_wall_flux_solar_rad", "west_exterior_wall_flux_solar_rad", "south_exterior_wall_flux_solar_rad"];
data, column_names = sumColumns(data, column_names, exterior_flux_names, "exterior_wall_flux_solar_rad");

exterior_flux_names = ["east_exterior_wall_flux_sky_rad", "north_exterior_wall_flux_sky_rad", "west_exterior_wall_flux_sky_rad", "south_exterior_wall_flux_sky_rad"];
data, column_names = sumColumns(data, column_names, exterior_flux_names, "exterior_wall_flux_sky_rad");

net_flux_names = ["east_interior_wall_flux", "north_interior_wall_flux", "west_interior_wall_flux", "south_interior_wall_flux", "floor_flux", "ceiling_flux"];
data, column_names = sumColumns(data, column_names, net_flux_names, "net_interior_flux");

net_flux_names = ["east_exterior_wall_flux", "north_exterior_wall_flux", "west_exterior_wall_flux", "south_exterior_wall_flux", "roof_flux"];
data, column_names = sumColumns(data, column_names, net_flux_names, "net_exterior_flux");

net_flux_names = ["east_exterior_wall_flux_cond", "north_exterior_wall_flux_cond", "west_exterior_wall_flux_cond", "south_exterior_wall_flux_cond", "roof_flux_cond"];
data, column_names = sumColumns(data, column_names, net_flux_names, "net_exterior_flux_cond");

net_flux_names = ["east_exterior_wall_flux_solar_rad", "north_exterior_wall_flux_solar_rad", "west_exterior_wall_flux_solar_rad", "south_exterior_wall_flux_solar_rad", "roof_flux_solar_rad"];
data, column_names = sumColumns(data, column_names, net_flux_names, "net_exterior_flux_solar_rad");

net_flux_names = ["east_exterior_wall_flux_sky_rad", "north_exterior_wall_flux_sky_rad", "west_exterior_wall_flux_sky_rad", "south_exterior_wall_flux_sky_rad", "roof_flux_sky_rad"];
data, column_names = sumColumns(data, column_names, net_flux_names, "net_exterior_flux_sky_rad");

net_flux_names = ["floor_flux_north_rad", "floor_flux_east_rad", "floor_flux_south_rad", "floor_flux_west_rad"];
data, column_names = sumColumns(data, column_names, net_flux_names, "net_floor_flux_rad");

net_flux_names = ["air_leakage", "air_ventilation"];
data, column_names = sumColumns(data, column_names, net_flux_names, "net_air_ventilation");

data, column_names = copyColumn(data, column_names, "interior_wall_flux", -1, "interior_wall_flux_neg")
print("new column names = ", column_names)
data, column_names = copyColumn(data, column_names, "floor_flux_cond", -1, "floor_flux_cond_neg");
data, column_names = copyColumn(data, column_names, "ceiling_flux", -1, "ceiling_flux_neg");
data, column_names = copyColumn(data, column_names, "interior_wall_flux", -1, "interior_wall_flux_neg");
data, column_names = copyColumn(data, column_names, "net_air_ventilation", -1, "net_air_ventilation_neg")

net_flux_names = ["interior_wall_flux_neg", "floor_flux_cond_neg", "ceiling_flux_neg", "interior_wall_flux_neg", "net_air_ventilation_neg", "window_conduction"];
data, column_names = sumColumns(data, column_names, net_flux_names, "net_interior_load");




# plot interior results
fig, axs = plt.subplots(3);
fig.set_size_inches(8, 6)
fig.subplots_adjust(left=0.15, bottom=0.1, right=0.9, top=0.9, wspace=0.4, hspace=0.6)

column_names1 = ["time", "interior_air_temp", "exterior_air_temp"]
column_names2 = ["time", "hvac_flux", "floor_flux", "ceiling_flux", "interior_wall_flux", "net_interior_flux"]
column_names3 = ["time", "interior_air_temp", "floor_temp", "interior_wall_temp", "ceiling_temp"]

#column_names3 = ["time", "hvac_flux", "roof_flux", "exterior_wall_flux"]


plotColumns(column_names, column_names1, data, "time (s)", "temperature (K)", fig, axs[0])
plotColumns(column_names, column_names2, data, "time (s)", "flux (W)", fig, axs[1])
plotColumns(column_names, column_names3, data, "time (s)", "temperature (K)", fig, axs[2])


makeLegend(fig, axs)
fig.suptitle("Simple House")
fig.savefig("simple_house_interior_results.png", dpi=600);

# plot interior flux breakdown
fig, axs = plt.subplots(3);
fig.set_size_inches(8, 6)
fig.subplots_adjust(left=0.15, bottom=0.1, right=0.9, top=0.9, wspace=0.4, hspace=0.6)

column_names1 = ["time", "floor_flux_cond", "interior_wall_flux", "ceiling_flux"]
column_names2 = ["time", "floor_flux_north_rad", "floor_flux_east_rad", "floor_flux_south_rad", "floor_flux_west_rad"]
column_names3 = ["time", "air_leakage", "air_ventilation", "window_conduction", "interior_loads"]

plotColumns(column_names, column_names1, data, "time (s)", "flux (W)", fig, axs[0])
plotColumns(column_names, column_names2, data, "time (s)", "flux (W)", fig, axs[1])
plotColumns(column_names, column_names3, data, "time (s)", "flux (W)", fig, axs[2])

makeLegend(fig, axs)
fig.suptitle("Simple House")
fig.savefig("simple_house_interior_fluxes.png", dpi=600);


# plot exterior
fig, axs = plt.subplots(3);
fig.set_size_inches(8, 6)
fig.subplots_adjust(left=0.15, bottom=0.1, right=0.9, top=0.9, wspace=0.4, hspace=0.6)

column_names1 = ["time", "interior_air_temp", "exterior_air_temp"]
column_names2 = ["time", "hvac_flux", "roof_flux", "exterior_wall_flux", "net_exterior_flux"]
column_names3 = ["time",  "exterior_wall_temp", "roof_temp"]

plotColumns(column_names, column_names1, data, "time (s)", "temperature (K)", fig, axs[0])
plotColumns(column_names, column_names2, data, "time (s)", "flux (W)", fig, axs[1])
plotColumns(column_names, column_names3, data, "time (s)", "temperature (K)", fig, axs[2])

makeLegend(fig, axs)
fig.suptitle("Simple House")
fig.savefig("simple_house_exterior_results.png", dpi=600);

# plot exterior flux breakdown
fig, axs = plt.subplots(3);
fig.set_size_inches(8, 6)
fig.subplots_adjust(left=0.15, bottom=0.1, right=0.9, top=0.9, wspace=0.4, hspace=0.6)

column_names1 = ["time", "roof_flux_cond", "exterior_wall_flux_cond", "net_exterior_flux_cond"]
column_names2 = ["time", "roof_flux_solar_rad", "exterior_wall_flux_solar_rad", "net_exterior_flux_solar_rad"]
column_names3 = ["time", "roof_flux_sky_rad", "exterior_wall_flux_sky_rad", "net_exterior_flux_sky_rad"]


plotColumns(column_names, column_names1, data, "time (s)", "flux (W)", fig, axs[0])
plotColumns(column_names, column_names2, data, "time (s)", "flux (W)", fig, axs[1])
plotColumns(column_names, column_names3, data, "time (s)", "flux (W)", fig, axs[2])

makeLegend(fig, axs)
fig.suptitle("Simple House")
fig.savefig("simple_house_exterior_fluxes.png", dpi=600);

# plot under ground temperatures
fig, axs = plt.subplots(2);
fig.set_size_inches(8, 6)
fig.subplots_adjust(left=0.15, bottom=0.1, right=0.9, top=0.9, wspace=0.4, hspace=0.6)

column_names1 = ["time", "found_bottom_temp", "found_insl_bottom_temp", "ground_bottom_temp"]
plotColumns(column_names, column_names1, data, "time (s)", "temperature (K)", fig, axs[0])

makeLegend(fig, axs)
fig.suptitle("Simple House")
fig.savefig("underground_temps.png", dpi=600);

# plot interior net flux
fig, axs = plt.subplots(2);
fig.set_size_inches(8, 6)
fig.subplots_adjust(left=0.15, bottom=0.1, right=0.9, top=0.9, wspace=0.4, hspace=0.6)

column_names1 = ["time", "hvac_flux", "interior_wall_flux_neg", "ceiling_flux_neg", "floor_flux_cond_neg", "net_air_ventilation_neg", "window_conduction", "net_interior_load"]
column_labels1 = ["time", "hvac", "wall", "ceiling", "floor", "air_exch", "window_cond", "net_load"]
plotColumnsWithLabels(column_names, column_names1, data, "time (s)", "flux (W)", fig, axs[0], column_labels1)

column_names1 = ["time", "interior_wall_flux_neg", "ceiling_flux_neg", "floor_flux_cond_neg", "net_air_ventilation_neg", "window_conduction"]
column_labels1 = ["time", "wall", "ceiling", "floor", "air_exch", "window_cond"]
plotColumnsWithLabels(column_names, column_names1, data, "time (s)", "flux (W)", fig, axs[1], column_labels1)

makeLegend(fig, axs)
fig.suptitle("Simple House")
fig.savefig("interior_net_flux.png", dpi=600);

# plot some things in Imperial units
fig, axs = plt.subplots(3);
fig.set_size_inches(8, 6)
fig.subplots_adjust(left=0.15, bottom=0.1, right=0.9, top=0.9, wspace=0.4, hspace=0.6)

data, column_names = convertTemperatures(data, column_names)
data, column_names = convertFluxesToTons(data, column_names)

column_names4 = ["time", "interior_air_temp_F"]
column_names5 = ["time", "exterior_air_temp_F"]
column_names6 = ["time", "hvac_flux_tons", "interior_wall_flux_tons", "exterior_wall_flux_tons"]
plotColumns(column_names, column_names4, data, "time (s)", "temperature (F)", fig, axs[0])
plotColumns(column_names, column_names5, data, "time (s)", "temperature (F)", fig, axs[1])
plotColumns(column_names, column_names6, data, "time (s)", "tons", fig, axs[2])

makeLegend(fig, axs)
fig.suptitle("Simple House")
fig.savefig("simple_house_imperial_results.png", dpi=600)

computeTotalHVACFluxes(column_names, data)

writeResults(column_names, data, fname + "_postprocessed");