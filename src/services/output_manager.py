class OutputManager:

    def output_line(self, line: str, is_title: bool = False, end_line="\n"):
        if is_title:
            line = self.generate_title(line)
        print(line, end=end_line)

    def generate_title(self, title: str):
        line_length = 50  # This arbitrary value can be changed. Move to ENV?

        # Calculate the padding length for centering
        title_padding = (line_length - len(title)) // 2

        # Pad the title line with spaces to make it centered
        title_line = title.center(line_length, '=')

        # Print the title and content
        return title_line

    def output_heading(self):
        self.output_line("", is_title=True)
        self.output_line("compAna: a tool for comparing annotations")
        self.output_line("", is_title=True)
        self.output_line("\n")


default_output_manager = OutputManager()
