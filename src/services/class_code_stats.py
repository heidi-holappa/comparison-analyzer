from services.output_manager import default_output_manager as output_manager


class ClassCodeStats:

    def __init__(self, gffcompare_db):
        self.gffcompare_db = gffcompare_db

    def compute_class_code_stats(self):
        output_manager.output_line("CLASS CODE STATISTICS", is_title=True)
        class_codes = {}
        for transcript in self.gffcompare_db.features_of_type('transcript'):
            if 'class_code' in transcript.attributes:
                class_code = transcript.attributes['class_code'][0]
                if not class_code in class_codes:
                    class_codes[class_code] = 0
                class_codes[class_code] += 1

        output_manager.output_line(
            f"{'class code':<18}| {'n':<10}")
        output_manager.output_line("", is_title=True, fill="-")
        for key, value in sorted(class_codes.items(), key=lambda item: item[1], reverse=True):
            output_manager.output_line(f"{key:<18}| {value:<10}")

        output_manager.output_line("", is_title=True)
