#!/usr/bin/env sh

# Run a simple test of the null libland library

# Determine the where this script lives
script_root=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd -P)

get_executable_name(){
  # The executable name for the script is based on the script name.  The script
  # name will follow the name t<number>-<name>.<extension> where <number is a
  # string of numbers, <name> is the name of the executable, and <extension>
  # is the extension of the scripting language.
  stringToConvert=${1#*t*-}
  stringToConvert=${stringToConvert%.*}
  echo $stringToConvert
}
# Helper variables for use in this script
dataDir=${script_root}/test_data
gridDir=${dataDir}/grids/CM2.1

# Create the run directory for this script, based on the name of this script.
# Delete the contents if the directory already exists.
runDir=${script_root}/$(basename ${0%.*})
if [ -e ${runDir} ]
then
  rm -rf ${runDir}
fi
mkdir -p ${runDir}
# Verify runDir exists, and is a directory
if [ ! -e ${runDir} -a ! -d ${runDir} ]
then
  echo "The run directory ($runDir) for $0 does not exist.  Exiting." 1>&2
  exit 1
fi
# Finish preparing the run directory
mkdir -p ${runDir}/INPUT
if [ ! -e ${runDir}/INPUT -a ! -d ${runDir} ]
then
  echo "Unable to setup run directory.  Exiting." 1>&2
  exit 1
fi

# Convert the grid files to netCDF files:
for f in ${gridDir}/*.cdl
do
  ncgen -o ${runDir}/INPUT/$(basename ${f%.cdl}.nc) ${f}
done

cat >>${runDir}/field_table<<EOF
  "TRACER", "land_mod", "sphum"
           "longname",     "specific humidity"
            "units",        "kg/kg" /
EOF

touch ${runDir}/input.nml

# Enter the run directory, and run the executable, and exit script with the
# exit code of the executable.
cd ${runDir}
../$(get_executable_name $0)
exit $?
