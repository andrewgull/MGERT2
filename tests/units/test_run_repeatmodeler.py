import os
import sys
import pytest
from unittest.mock import MagicMock, patch, mock_open

# Add the project root to sys.path to allow importing from workflow.scripts
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from workflow.scripts.run_repeatmodeler import run_repeatmodeler_logic


def test_run_repeatmodeler_logic_success():
    """Test successful execution of RepeatModeler logic."""
    database = "test_db"
    threads = 4
    output_file = "output/families.fa"
    log_file = "logs/repeatmodeler.log"

    with patch(
        "workflow.scripts.run_repeatmodeler.subprocess.Popen"
    ) as mock_popen, patch(
        "workflow.scripts.run_repeatmodeler.glob.glob"
    ) as mock_glob, patch(
        "workflow.scripts.run_repeatmodeler.os.path.getmtime"
    ) as mock_getmtime, patch(
        "workflow.scripts.run_repeatmodeler.os.path.exists"
    ) as mock_exists, patch(
        "workflow.scripts.run_repeatmodeler.os.makedirs"
    ) as mock_makedirs, patch(
        "workflow.scripts.run_repeatmodeler.shutil.copy2"
    ) as mock_copy, patch(
        "workflow.scripts.run_repeatmodeler.shutil.rmtree"
    ) as mock_rmtree, patch(
        "workflow.scripts.run_repeatmodeler.logging.basicConfig"
    ), patch(
        "builtins.open", mock_open()
    ):
        # Mock subprocess behavior
        mock_process = MagicMock()
        mock_process.returncode = 0
        mock_popen.return_value = mock_process

        # Mock glob and directory finding
        mock_glob.return_value = ["RM_1", "RM_2"]
        # RM_2 is newer
        mock_getmtime.side_effect = lambda x: 100 if x == "RM_1" else 200
        mock_exists.return_value = True

        run_repeatmodeler_logic(database, threads, output_file, log_file)

        # Assertions
        mock_popen.assert_called_once()
        command = mock_popen.call_args[0][0]
        assert "RepeatModeler" in command
        assert "-threads" in command
        assert "4" in command
        assert "-database" in command
        assert "test_db" in command

        mock_process.wait.assert_called_once()
        mock_glob.assert_called_with("RM_*")
        mock_copy.assert_called_once_with("RM_2/consensi.fa.classified", output_file)
        mock_rmtree.assert_called_once_with("RM_2")


def test_run_repeatmodeler_logic_failure():
    """Test RepeatModeler failure handling (non-zero exit code)."""
    database = "test_db"
    threads = 4
    output_file = "output/families.fa"
    log_file = "logs/repeatmodeler.log"

    with patch(
        "workflow.scripts.run_repeatmodeler.subprocess.Popen"
    ) as mock_popen, patch(
        "workflow.scripts.run_repeatmodeler.logging.basicConfig"
    ), patch(
        "builtins.open", mock_open()
    ), patch(
        "workflow.scripts.run_repeatmodeler.sys.exit"
    ) as mock_exit:
        mock_process = MagicMock()
        mock_process.returncode = 1
        mock_popen.return_value = mock_process
        mock_exit.side_effect = SystemExit

        with pytest.raises(SystemExit):
            run_repeatmodeler_logic(database, threads, output_file, log_file)

        mock_exit.assert_called_once_with(1)


def test_run_repeatmodeler_no_dir():
    """Test handling when no RM_ directory is found."""
    database = "test_db"
    threads = 4
    output_file = "output/families.fa"
    log_file = "logs/repeatmodeler.log"

    with patch(
        "workflow.scripts.run_repeatmodeler.subprocess.Popen"
    ) as mock_popen, patch(
        "workflow.scripts.run_repeatmodeler.logging.basicConfig"
    ), patch(
        "builtins.open", mock_open()
    ), patch(
        "workflow.scripts.run_repeatmodeler.glob.glob"
    ) as mock_glob, patch(
        "workflow.scripts.run_repeatmodeler.sys.exit"
    ) as mock_exit:
        mock_process = MagicMock()
        mock_process.returncode = 0
        mock_popen.return_value = mock_process
        mock_glob.return_value = []
        mock_exit.side_effect = SystemExit

        with pytest.raises(SystemExit):
            run_repeatmodeler_logic(database, threads, output_file, log_file)

        mock_exit.assert_called_once_with(1)


def test_run_repeatmodeler_missing_output_file():
    """Test handling when the expected output file is missing in the RM_ directory."""
    database = "test_db"
    threads = 4
    output_file = "output/families.fa"
    log_file = "logs/repeatmodeler.log"

    with patch(
        "workflow.scripts.run_repeatmodeler.subprocess.Popen"
    ) as mock_popen, patch(
        "workflow.scripts.run_repeatmodeler.logging.basicConfig"
    ), patch(
        "builtins.open", mock_open()
    ), patch(
        "workflow.scripts.run_repeatmodeler.glob.glob"
    ) as mock_glob, patch(
        "workflow.scripts.run_repeatmodeler.os.path.getmtime"
    ) as mock_getmtime, patch(
        "workflow.scripts.run_repeatmodeler.os.path.exists"
    ) as mock_exists, patch(
        "workflow.scripts.run_repeatmodeler.sys.exit"
    ) as mock_exit:
        mock_process = MagicMock()
        mock_process.returncode = 0
        mock_popen.return_value = mock_process
        mock_glob.return_value = ["RM_1"]
        mock_getmtime.return_value = 100
        mock_exists.return_value = False  # Output file missing
        mock_exit.side_effect = SystemExit

        with pytest.raises(SystemExit):
            run_repeatmodeler_logic(database, threads, output_file, log_file)

        mock_exit.assert_called_once_with(1)
