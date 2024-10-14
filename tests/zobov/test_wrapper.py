from unittest import mock

from voidfindertk.zobov import _wrapper


# @mock.patch("sh.Command")  # Patching sh.Command
def test_vozinit(load_mock_data):
    mock_data = load_mock_data("zobov", "vozinit.jpkl")

    # extract the output value and left alone the inputs
    expected_output = mock_data.pop("output")

    # put "." as working dir
    mock_data["work_dir_path"] = "."

    # mock the instance of command mock
    mock_vozinit = mock.Mock(return_value=expected_output)

    # mock the class and run
    with mock.patch("sh.Command", return_value=mock_vozinit):
        output = _wrapper.run_vozinit(**mock_data)

    assert output == expected_output
    mock_vozinit.assert_called_once_with(
        str(mock_data["input_file_path"]),
        str(mock_data["buffer_size"]),
        str(mock_data["box_size"]),
        str(mock_data["number_of_divisions"]),
        str(mock_data["executable_name"]),
    )
