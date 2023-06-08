from datetime import datetime


class OutputManager:

    def output_line(self,
                    line: str,
                    is_title: bool = False,
                    end_line: str = "\n",
                    fill: str = '=',
                    additional_line_breaks: int = 0,
                    is_info: bool = False,
                    is_error: bool = False,
                    title_line_length: int = 50):
        if is_title:
            line = self.generate_title(line, fill, title_line_length)
        if is_info:
            line = f"INFO: [{datetime.now().strftime('%Y-%m-%d %H:%M')}]: {line}"
        if is_error:
            line = f"ERROR: [{datetime.now().strftime('%Y-%m-%d %H:%M')}]: {line}"
        print(line, end=end_line)
        if additional_line_breaks:
            print("\n" * additional_line_breaks, end="")

    def generate_title(self, title: str, fill: str, line_length: int):
        if len(title):
            title = " " + title + " "
        title_line = (title).center(line_length, fill)
        return title_line

    def output_heading(self):
        self.output_line("", is_title=True)
        self.output_line(
            "compAna: a tool for comparing annotations", fill=' ', is_title=True)
        self.output_line("", is_title=True)

    def output_footer(self):
        self.output_line("Pipeline finished.", is_info=True)
        self.output_line("END", is_title=True)


default_output_manager = OutputManager()
