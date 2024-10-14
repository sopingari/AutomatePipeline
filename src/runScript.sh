
SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

source ${SCRIPT_DIR}/miniconda3/bin/activate base

export exit_code=0
python -m cc3d.run_script $* --current-dir=${SCRIPT_DIR}
exit_code=$?

cd ${SCRIPT_DIR}
exit ${exit_code}