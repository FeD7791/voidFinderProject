import sh 

from ..utils import chdir

def _run_main(*,main_path,vf_param_path,work_dir_path):
    main = sh.Command('main.x',search_paths=[main_path])
    args = [str(vf_param_path)]
    with chdir(work_dir_path):
        output = main(*args)
    return output

    
