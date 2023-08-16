import unittest
import pytest
from datetime import datetime
from services.output_manager import default_output_manager as output_manager


class TestOutputManager(unittest.TestCase):

    @pytest.fixture(autouse=True)
    def capsys(self, capsys):
        self.capsys = capsys

    def test_an_empty_string_title_generates_a_line_of_selected_character(self):
        output_manager.output_line({
            "line": "",
            "is_title": True,
            "title_line_length": 10
        })
        captured = self.capsys.readouterr()
        assert captured.out == "==========\n"

    def test_a_string_title_generates_a_line_of_selected_character(self):
        output_manager.output_line({
            "line": "test",
            "is_title": True,
            "title_line_length": 10
        })
        captured = self.capsys.readouterr()
        assert captured.out == "== test ==\n"

    def test_a_string_title_generates_a_line_of_selected_character_with_one_additional_line_break(self):
        output_manager.output_line({
            "line": "test",
            "is_title": True,
            "title_line_length": 10,
            "additional_line_breaks": 1
        })
        captured = self.capsys.readouterr()
        assert captured.out == "== test ==\n\n"

    def test_a_string_title_generates_a_line_of_selected_character_with_two_additional_line_breaks(self):
        output_manager.output_line({
            "line": "test",
            "is_title": True,
            "title_line_length": 10,
            "additional_line_breaks": 2
        })
        captured = self.capsys.readouterr()
        assert captured.out == "== test ==\n\n\n"

    def test_a_string_with_empty_end_line_generates_a_line_with_no_end_line(self):
        output_manager.output_line({
            "line": "test",
            "end_line": ""
        })
        captured = self.capsys.readouterr()
        assert captured.out == "test"

    def test_a_string_w_arg_info_generates_a_line_with_info_prefix(self):
        output_manager.output_line({
            "line": "test",
            "is_info": True
        })
        captured = self.capsys.readouterr()
        assert "INFO:" in captured.out

    def test_a_string_w_arg_error_generates_a_line_with_error_prefix(self):
        output_manager.output_line({
            "line": "test",
            "is_error": True
        })
        captured = self.capsys.readouterr()
        assert "ERROR:" in captured.out

    def test_a_string_w_arg_is_info_generates_a_line_with_info_prefix_and_date(self):
        output_manager.output_line({
            "line": "test",
            "is_info": True
        })
        captured = self.capsys.readouterr()
        assert "INFO: [{}]:".format(
            datetime.now().strftime('%Y-%m-%d %H:%M')) in captured.out

    def test_generating_heading_creates_lines_of_string(self):
        output_manager.output_heading()
        captured = self.capsys.readouterr()
        assert "compAna" in captured.out

    def test_generating_footer_creaters_lines_of_string(self):
        run_time_str = "test"
        output_manager.output_footer(run_time_str)
        captured = self.capsys.readouterr()
        assert "compAna" in captured.out
