def compute_class_code_stats(gffcompare_db):
    print("==========CLASS CODE STATISTICS==========")
    class_codes = {}
    print('\nComputing statistics for class codes...\n')
    for transcript in gffcompare_db.features_of_type('transcript'):
        if 'class_code' in transcript.attributes:
            class_code = transcript.attributes['class_code'][0]
            if not class_code in class_codes:
                class_codes[class_code] = 0
            class_codes[class_code] += 1

    print(f"{'class code':<18}| {'n':<10}")
    print('-' * 30)
    for key, value in sorted(class_codes.items(), key=lambda item: item[1], reverse=True):
        print(f"{key:<18}| {value:<10}")

    print("\n=========================================")
