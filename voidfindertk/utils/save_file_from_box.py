import pandas as pd


def box_to_csv(box, file_path, **kwargs):
    kwargs.setdefault("path_or_buf", file_path)
    kwargs.setdefault("delim_whitespace", True)
    kwargs.setdefault("header", False)
    kwargs.setdefault("index", False)
    df = pd.DataFrame(box.__dict__)
    df.to_csv(**kwargs)
