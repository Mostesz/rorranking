import argparse
import os
from utils import get_all_files_paths


def merge_files(output_path, data_path):
    files_paths = get_all_files_paths(data_path)
    for file_path in files_paths:
        output_file_path = os.path.join(output_path, file_path[file_path.index('/') + 1:])
        if not os.path.exists(os.path.dirname(output_file_path)):
            os.makedirs(os.path.dirname(output_file_path))
        with open(output_file_path, "a") as output_file:
            input_file_path = os.path.join(data_path, file_path)
            lines_number = 0
            with open(input_file_path, "r") as input_file:
                for line in input_file:
                    output_file.write(line)
                    lines_number += 1


def verify_files(path):
    files_paths = get_all_files_paths(path)
    for file_path in files_paths:
        with open(os.path.join(path, file_path), "r") as input_file:
            lines_number = 0
            for line in input_file:
                if len(line.strip()) > 0:
                    lines_number += 1
            if lines_number < 100:
                print "Error for path: %s" % file_path


def merge_data(output_path, data_path):
    merge_files(output_path, data_path)
    verify_files(output_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("output_path", help="output path")
    parser.add_argument("data_path", help="path to results of experiments")
    args = parser.parse_args()

    merge_data(args.output_path, args.data_path)
