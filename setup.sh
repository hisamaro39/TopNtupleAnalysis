SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  CURRDIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$CURRDIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
CURRDIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CURRDIR/../TuDoBase/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CURRDIR/../TuDoBase/Root