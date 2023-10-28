
import pandas as pd

from . import box


def read_table(path_or_buffer, **kwargs):

    kwargs.setdefault("sep", r"\s+")
    kwargs.setdefault("usecols", [0, 1, 2, 3, 4, 5, 6])
    data = pd.read_csv(path_or_buffer, **kwargs)
    col_number = len(data.columns)

    ##Check for 7 columns
    if col_number != 7:
        raise ValueError(
            "There are not enough columns to create the coordinates of a box. "
            f"Found {col_number} expected 7"
            )
    ##Check for null values
    check_values = data.notnull().values.all()
    if not check_values:
        raise TypeError(f'There are:' 
                        f'{data.isnull().sum().sum()} null or missing values')

    the_box = box.Box(
        x=data.iloc[:, 0],
        y=data.iloc[:, 1],
        z=data.iloc[:, 2],
        vx=data.iloc[:, 3],
        vy=data.iloc[:, 4],
        vz=data.iloc[:, 5],
        m =data.iloc[:, 6],
    )
    return the_box