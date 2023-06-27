import os
from datetime import datetime
from config import LOG_FILE_DIR


class OutputManager:

    def __init__(self):
        self.log_output = []

    def output_line(self, config: dict):
        """Outputs given line to console. Can be used to output titles, info, errors and more. A simple custom logger.

        Args:
            config (dict): Configuration options for the output line.
        """
        line = config.get('line', '')
        is_title = config.get('is_title', False)
        end_line = config.get('end_line', '\n')
        fill = config.get('fill', '=')
        additional_line_breaks = config.get('additional_line_breaks', 0)
        is_info = config.get('is_info', False)
        is_error = config.get('is_error', False)
        title_line_length = config.get('title_line_length', 50)
        save_to_log = config.get('save_to_log', True)
        if is_title:
            line = self.generate_title(line, fill, title_line_length)
        if is_info:
            line = f"INFO: [{datetime.now().strftime('%Y-%m-%d %H:%M')}]: {line}"
        if is_error:
            line = f"ERROR: [{datetime.now().strftime('%Y-%m-%d %H:%M')}]: {line}"
        print(line, end=end_line)
        if additional_line_breaks:
            print("\n" * additional_line_breaks, end="")
        if save_to_log:
            self.log_output.append(line + "\n")

    def generate_title(self, title: str, fill: str, line_length: int):
        if len(title):
            title = " " + title + " "
        title_line = (title).center(line_length, fill)
        return title_line

    def output_heading(self):
        self.output_line({
            "line": "",
            "is_title": True
        })
        self.output_line(
            {
                'line': "compAna: a tool for comparing annotations",
                'fill': ' ',
                'is_title': True
            }
        )
        self.output_line(
            {
                "line": "",
                "is_title": True
            }
        )

    def output_footer(self):
        self.output_line({
            "line": "",
            "is_title": True
        })
        self.output_line({
            "line": "Pipeline finished.",
            "is_info": True
        })
        self.output_line({
            "line": "Thank you for using compAna",
            "is_info": True
        })

    def write_log_file(self):
        if not self.log_output:
            return
        file_path = os.path.join(LOG_FILE_DIR, 'stdout.log')
        with open(file_path, 'a', encoding='utf-8') as file:
            file.writelines(self.log_output)
        self.output_line({
            "line": "stdout appended to " + file_path,
            "is_info": True
        })


default_output_manager = OutputManager()
