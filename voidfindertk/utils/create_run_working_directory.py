import pathlib

import datetime as dt

def create_run_work_dir(workdir_path):
        """
        This method will create a temporal directory inside workdir
        """
        timestamp = dt.datetime.now(dt.timezone.utc).isoformat()
        run_work_dir = pathlib.Path(
            tempfile.mkdtemp(suffix=timestamp, dir=work_dir_path)
        )
        return run_work_dir