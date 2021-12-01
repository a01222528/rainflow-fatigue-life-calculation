# %%
import pandas
import numpy

filename = "Test6"

def compare_plots(base,modified,params):
    plot = pandas.merge(
        base,
        modified,
        how="outer",
        left_index=True,
        right_index=True
        ).interpolate(method="linear").plot(**params,grid=True,zorder=2)
    plot.get_figure().savefig("rf_"+modified.keys()[0]+".png")

def get_number_of_decimals( flt ):
    for i in range(10):
        if flt*(10**i) == round(flt*(10**i)):
            break
    else:
        raise Exception(flt, "is too long to get number of decimals.")
    return i

plot_parameters = {
    "xlabel" : "timestamp",
    "ylabel" : "MPa",
    "xlim" : (320,345),
    "ylim" : (1,8)
}

algorithm_parameters = {
    "hysteresis_gate"   : 0.1,
    "binning" : {
        "value"         : 0.5,
        "decimals"      : None
    }
}
algorithm_parameters["binning"]["decimals"] = get_number_of_decimals( algorithm_parameters["binning"]["value"] )


# %%
# Data Preprocessing
df = pandas.read_csv(filename+".csv")[7:]

df.columns="ts,time,stress".split(",")
df[["time","stress"]] = df[["time","stress"]].apply(pandas.to_numeric)

df["stress"] = df["stress"]*205000/1000000

stress = df["stress"].values.tolist()
stress_serie = pandas.DataFrame(stress,columns=["stress"])

total_time = df["time"].max()-df["time"].min()

print("Min stress:", stress_serie.min()[0])
print("Max stress:", stress_serie.max()[0])
print("Measurements:", len(stress_serie))
print("Total time (sec):", total_time)
print("Range:", stress_serie.max()[0]-stress_serie.min()[0])
print("Maximum bin value:", (stress_serie.max()[0]-stress_serie.min()[0])/64)

stress_serie.plot( title = "Complete Series", **dict(item for item in plot_parameters.items() if item[0] in ["xlabel","ylabel"]) )
stress_serie[320:345].plot( title = "Head Series", **dict(item for item in plot_parameters.items() if item[0] in ["xlabel","ylabel"]) )

# %%
# Hysteresis Filtering
def hysteresis_filtering( data, gate ):
    ftr = data.copy(deep=True)
    i = 0
    to_drop = list()
    while i < len(data):
        j = 1
        try:
            while abs(ftr["stress"][i+j] - ftr["stress"][i]) < gate:
                to_drop.append(i+j)
                j += 1
            i = i+j
        except:
            break
    ftr.drop(to_drop,inplace=True)

    print("Input length: {}\nOutput length: {}\nRemoved: {}".format(
        len(data),
        len(ftr),
        len(data)-len(ftr)
    ))

    return ftr

stress_hysteresis_filtered = hysteresis_filtering( stress_serie, algorithm_parameters["hysteresis_gate"] )
compare_plots(stress_serie,stress_hysteresis_filtered.rename(columns = {"stress":"hysteresis_filtered"}),{"title":"Hysteresis Filtering",**plot_parameters})

# %%
# Peak-Valley Filtering
def peak_valley_filtering( data ):
    ftr = data.copy(deep=True)
    i = 1
    idx = data["stress"].keys()
    to_drop = list()
    orientation_prev = "down" if ftr["stress"][idx[1]] < ftr["stress"][idx[0]] else "up"
    while i < len(data)-1:
        orientation = "down" if ftr["stress"][idx[i+1]] < ftr["stress"][idx[i]] else "up"
        if orientation_prev == orientation:
            to_drop.append(idx[i])
        else:
            orientation_prev = orientation
        i=i+1
  
    ftr.drop(to_drop,inplace=True)

    print("Input length: {}\nOutput length: {}\nRemoved: {}".format(
        len(data),
        len(ftr),
        len(data)-len(ftr)
    ))

    return ftr

stress_peak_valley_filtering = peak_valley_filtering( stress_hysteresis_filtered )
compare_plots(stress_serie,stress_peak_valley_filtering.rename(columns = {"stress":"peak_valley_filtered"}),{"title":"Peak-Valley Filtering",**plot_parameters})

# %%
# Binning
def binning( data, bin ):
    ftr = data.copy(deep=True)
    i = 0
    idx = data["stress"].keys()
    while i < len(data):
        ftr["stress"][idx[i]] = (ftr["stress"][idx[i]]//bin["value"])*bin["value"]+bin["value"]/2
        i=i+1

    return ftr

stress_bin = binning( stress_peak_valley_filtering, algorithm_parameters["binning"] )
compare_plots(stress_serie,stress_bin.rename(columns = {"stress":"binned"}),{"title":"Binning",**plot_parameters})

input_len = len(set(i[0] for i in stress_peak_valley_filtering.values))
output_len = len(set([i[0] for i in stress_bin.values]))

print("Input length: {}\nOutput length: {}\nRemoved: {}".format(
    input_len,
    output_len,
    input_len-output_len
    ))

# %%
# Four Point Counting Method
def four_point_counting_method( data ):

    matrix_len = numpy.arange(data.min()[0],data.max()[0]+algorithm_parameters["binning"]["value"],algorithm_parameters["binning"]["value"]).round(5)
    print("Matrix size: ", len(matrix_len))

    matrix_dict = dict()
    for f in matrix_len:
        matrix_dict[f] = dict()
        for t in matrix_len:
            matrix_dict[f][t] = 0

    ftr = data.copy(deep=True)
    idx = data["stress"].keys()

    s = [0, 1, 2, 3]

    while s[3] < len(idx):                
        Z = [ ftr["stress"][idx[i]].round(5) for i in s ]

        if  (Z[0] <= Z[1] and Z[0] <= Z[2] and Z[1] <= Z[3] and Z[2] <= Z[3]) or \
            (Z[3] <= Z[1] and Z[3] <= Z[2] and Z[1] <= Z[0] and Z[2] <= Z[0]):
            
            # Complete Cycle
            matrix_dict[Z[1]][Z[2]] += 1

            idx = idx.drop(idx[s[1]])
            idx = idx.drop(idx[s[2]])

            s = [0, 1, 2, 3]

        else:
            # Incomplete Cycle
            s[0] = s[1]
            s[1] = s[2]
            s[2] = s[3]
            s[3] = s[3] + 1

    ftr.drop(set(data["stress"].keys())-set(idx),inplace=True)

    return matrix_dict,ftr

matrix_dict,residue = four_point_counting_method( stress_bin )

residue_plot_parameters = {
    "title":"Residue",
    "xlabel" : "timestamp",
    "ylabel" : "MPa"
}
compare_plots(stress_serie,residue.rename(columns = {"stress":"residue"}),residue_plot_parameters)
print("Number of residue elements",len(residue))

# %%
def plot_rainflow_matrix(matrix_dict, title, outname):
    lims = (-30,30)
    rainflow_matrix = list()
    for f in matrix_dict:
        for t in matrix_dict[f]:
            rainflow_matrix.append( [f,t,matrix_dict[f][t]] )

    last_df = pandas.DataFrame(rainflow_matrix, columns=["from","to","count"])
    plot = last_df.plot.scatter(title=title,x="to",y="from",c="count",xlabel="to (MPa)",ylabel="from (MPa)",grid=True,colormap="YlOrRd",sharex=False,xlim=lims,ylim=lims,zorder=2)
    plot.get_figure().savefig(filename+"/"+outname)
    return rainflow_matrix

rainflow_matrix = plot_rainflow_matrix( matrix_dict, "Rainflow Matrix", "rf_rainflow_matrix.png")
print("Rainflow matrix size:", len(rainflow_matrix))

# %%
# Rainflow information
rainflow_info_table = list()
for f,t,c in rainflow_matrix:
    rainflow_info_table.append([c, (f+t)/2, abs(t-f) ])
rainflow_info_table = numpy.array([info for info in rainflow_info_table if info[0] != 0])

n, S_m, S_a = rainflow_info_table[:,0],rainflow_info_table[:,1],rainflow_info_table[:,2]

# %%
# Stress-life methodology
import numpy as np
import math

# Material Properties
SRI1 = 828.674176498713 
b1 = -0.0760625751181736 
b2 = -0.0395348487609183 
Nc1 = 5e8 
UTS = 350 
El = 0 

S = (S_a)/(1-S_m/UTS)
sorted_counts_and_stress = sorted([[i,s] for i,s in zip(n,S)],key=lambda x: x[1])

a_df = pandas.DataFrame(sorted_counts_and_stress)
a_array = a_df[a_df[1] != 0].values

n,S = a_array[:,0],a_array[:,1]

kf = 1 
N1e6 = 1e6**b1*SRI1/2 
b1 = (math.log10(N1e6/kf)-math.log10(SRI1/2))/(math.log10(1e6)-math.log10(1)) 

k = (1/(SRI1/2))**(1/b1) 
m = 1/b1 

N = k*(S)**m 

Nm1e3 = len(N[np.where(N<1000)[0]]) 
Nm1e3 = len(N) - Nm1e3 + 1 

SRI10 = 2*UTS 
N1e3 = 1e3**b1*SRI1/2 
b0 = (math.log10(N1e3)-math.log10(SRI10/2))/(math.log10(1e3)-math.log10(1)) 
m00 = 1/b0 
k0 = (1/(SRI10/2))**(1/b0) 

eS_UTS = len(S[np.where(S<math.ceil(UTS))[0]])+1

N[Nm1e3:eS_UTS] = k0*(S[Nm1e3:eS_UTS])**m00 
N[eS_UTS:] = 1 

if(b2 != 0):
    m22 = 1/b2 
    SNc1 = (Nc1/k)**(1/m) 
    k2 = Nc1/SNc1**m22 
    NmNc1 = len(N[np.where(N>Nc1)[0]]) 
    N[:NmNc1] = k2*(S[:NmNc1])**m22 

# %%
# Calculate Damage
total_damage = n/N
print("Total damage (%): ", total_damage.sum())

# %%
# Expected Life
damage_per_second = i_damage.sum() / total_time 
expected_life_seconds = 1 / damage_per_second

print("Expected life in seconds: {}".format(expected_life_seconds.round(2)))
print("Expected life in hours:   {}".format((expected_life_seconds/3600).round(2)))
print("Expected life in years:   {}".format((expected_life_seconds/31536000).round(5)))