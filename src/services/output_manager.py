class OutputManager:

    def output_line(self, line: str, is_title: bool = False, end_line="\n", fill='=', additional_line_breaks=0):
        if is_title:
            line = self.generate_title(line, fill)
        print(line, end=end_line)
        if additional_line_breaks:
            print("\n" * additional_line_breaks)

    def generate_title(self, title: str, fill: str):
        line_length = 50  # This arbitrary value can be changed. Move to ENV?
        title_line = (" " + title + " ").center(line_length, fill)
        return title_line

    def output_heading(self):
        self.output_line("", is_title=True)
        self.output_line(
            "compAna: a tool for comparing annotations", fill=' ', is_title=True)
        self.output_line("", is_title=True, additional_line_breaks=2)


default_output_manager = OutputManager()
