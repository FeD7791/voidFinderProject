from voidfindertk.zobov import _wrapper
from unittest import TestCase, mock
from unittest.mock import MagicMock
import os
import pathlib
from io import StringIO
import numpy as np





@mock.patch('sh.Command')  # Patching sh.Command
@mock.patch('os.chdir',return_value=MagicMock)       # Mock chdir context manager
def test_run_vozinit(mock_chdir, mock_command):
    # Mock the vozinit command object
    mock_vozinit = MagicMock()
    mock_vozinit.return_value = "mocked output"
    
    # Mock what sh.Command returns
    mock_command.return_value = mock_vozinit

    # Define dummy input parameters
    vozinit_dir_path = "/fake/path"
    input_file_path = "input.txt"
    buffer_size = 1024
    box_size = 50
    number_of_divisions = 10
    executable_name = "vozinit_exec"
    work_dir_path = "/fake/workdir"

    # Run the function under test
    output = _wrapper.run_vozinit(
        vozinit_dir_path=vozinit_dir_path,
        input_file_path=input_file_path,
        buffer_size=buffer_size,
        box_size=box_size,
        number_of_divisions=number_of_divisions,
        executable_name=executable_name,
        work_dir_path=work_dir_path,
    )

    # Assertions
    mock_command.assert_called_once_with("vozinit", search_paths=[vozinit_dir_path])
    mock_vozinit.assert_called_once_with(
        input_file_path, str(buffer_size), str(box_size), str(number_of_divisions), executable_name
    )
    assert output == "mocked output"


