import os


def get_all_files_paths(base_directory):
    paths = []
    for root, directories, file_names in os.walk(base_directory):
        for file_name in file_names:
            paths.append(os.path.join(root, file_name)[len(base_directory.rstrip('/'))+1:])
    return paths
