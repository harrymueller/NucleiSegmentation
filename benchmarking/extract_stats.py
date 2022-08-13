
# libs
import pandas as pd
import numpy as np
import datetime
import sys

####################################
# FUNCTIONS: used for .apply(...)
####################################
# function to return ram in GB
def format_rss(x):
	if (x[-1] == "G"):
		return np.double(x[:-1])
	else:
		return np.double(x[:-1])/1000

# function to convert timedelta to nicely formatted string
def timedelta_to_str(x):
	seconds = int(x / np.timedelta64(1, 's'))
	return str(datetime.timedelta(seconds = round(seconds)))



####################################
# PIDS
####################################
def get_pids(pid_path):
	pids = pd.read_table(pid_path, sep=" ")

	# removing rows w/ names
	pids = pids.loc[pids["PID"] != "PID"].reset_index()

	# adding new columns to pid for statistics calculated later
	empty = np.empty((len(pids), 6))
	empty[:] = np.nan

	pids[["num_cells", "Time", "cpu_max", "cpu_avg", "ram_max", "ram_avg"]] = empty

	pids["PID"] = pids["PID"].astype(np.int32)

	for sample_id in NUM_CELLS:
		pids.loc[pids["SAMPLE_ID"] == sample_id, "num_cells"] = NUM_CELLS[sample_id]

	pids["num_cells"] = pids["num_cells"].astype(np.int32)
	pids = pids[pids.columns[1:]]

	return pids

####################################
# USAGE
####################################
def get_usage(usage_path):
	usage = [["Time", "UID", "PID", "%usr", "%system", "%guest", "%wait", "%CPU", "CPU", "minflt/s", "majflt/s", "VSZ", "RSS", "%MEM", "Command"]]

	# reading in from file
	with open(usage_path, "r") as f:
		f.readline()
		lines = f.readlines()

		for l in lines:
			if (l != "\n" and l[0] != "#" and l[0] != "L"):
				usage.append(l.replace("\n", "").replace(" AM", "AM").replace(" PM", "PM").split())
		
	# converting to pandas df and selecting only important columns
	usage = pd.DataFrame(usage[1:], columns=usage[0])
	usage = usage[["Time", "PID", "%CPU", "RSS"]]

	# casting + pd.Timedelta(1, unit='d')
	usage["Time"] = [pd.Timedelta(x) for x in usage["Time"]]
	usage["PID"] = usage["PID"].astype(np.int32)
	usage["%CPU"] = usage["%CPU"].apply(lambda x: x[:-1]).astype(np.double)
	usage["RSS"] = usage["RSS"].apply(lambda x: format_rss(x)).astype(np.double)

	return usage

####################################
# adds a day to the time unit if a time is less than the previous time
####################################
def add_days(usage):
	for i in range(1, len(usage["Time"])):
		if (usage["Time"][i] < usage["Time"][i-1]): # if new day
			for j in range(i, len(usage["Time"])):
				usage.at[j, "Time"] = usage.at[j, "Time"] + pd.Timedelta(1, unit="d")

	return usage

####################################
# getting stats from usage
####################################
def stats(usage):
	return {
		"time": usage.iloc[-1]["Time"] - usage.iloc[0]["Time"],
		"cpu_max": max(usage["%CPU"]),
		"cpu_avg": np.mean(usage["%CPU"]),
		"ram_max": max(usage["RSS"]),
		"ram_avg": np.mean(usage["RSS"])
	}

####################################
# print usage
####################################
def print_stats(cmd, usage_stats, output):
	out = "Time: \t\t{}\nMax CPU: \t{}%\nAvg. CPU: \t{}%\nMax RAM: \t{}GB\nAvg. RAM: \t{}GB".format(
		str(usage_stats["time"]).split(" ")[-1], # may have problems if processes > day
		usage_stats["cpu_max"],
		round(usage_stats["cpu_avg"], 1),
		round(usage_stats["ram_max"], 1),
		round(usage_stats["ram_avg"], 1)
	)
	out = "> {}\n{}".format(cmd,out)
	with open(output, "w") as f:
		f.write(out)

	print(out)
	print("\n")

####################################
# saving df to file
####################################
def save_df_to_file(tallies, pids, file):
	# converting timedelta to readable string
	pids["Time"] = pids["Time"].apply(lambda x: timedelta_to_str(x)).astype(str)


	with pd.ExcelWriter(file) as writer:
		tallies.to_excel(writer, sheet_name="Method Tallies", index=False)  
		pids.to_excel(writer, sheet_name="Methods and Sample IDs", index=False)  



####################################
# MAIN
####################################
def main():
	path = sys.argv[-2]
	output = sys.argv[-1]
	cmd = " ".join(sys.argv[1:-2])

	# input data
	usage = get_usage(path)
	usage = add_days(usage)

	usage_stats = stats(usage)
	print_stats(cmd, usage_stats, output)

if __name__ == "__main__":
	main()