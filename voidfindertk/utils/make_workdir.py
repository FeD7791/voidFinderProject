import datetime as dt
import pathlib
import tempfile


def create_run_work_dir(*, workdir_path):
    """
    This method will create a temporal directory inside the working
    directory of the ZobovVF class workdir.

    Returns
    -------
        run_work_dir: pathlib.Path
            path of the work directoty
    """
    timestamp = dt.datetime.now(dt.timezone.utc).isoformat()
    run_work_dir = pathlib.Path(
        tempfile.mkdtemp(suffix=timestamp, dir=workdir_path)
    )
    return run_work_dir
