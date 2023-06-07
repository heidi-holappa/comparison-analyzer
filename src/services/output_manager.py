class OutputManager:

    def output_line(self, line: str, is_title: bool = False, end_line="\n"):
        if is_title:
            line = self.generate_title(line)
        print(line, end=end_line)

    def generate_title(self, title: str):
        line_length = 40  # This arbitrary value can be changed. Move to ENV?

        # Calculate the padding length for centering
        title_padding = (line_length - len(title)) // 2

        # Pad the title line with spaces to make it centered
        title_line = title.center(line_length, '=')

        # Print the title and content
        return title_line


default_output_manager = OutputManager()
