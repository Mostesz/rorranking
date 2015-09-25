import os
import shutil
import unittest
from hamcrest import *
from os import path
from merge_results import merge_files


class MergeResultsTest(unittest.TestCase):
    def test_files_merging(self):
        input_path = path.join(path.dirname(__file__), 'fixtures', 'merge_results', 'input')
        output_path = path.join(path.dirname(__file__), 'fixtures', 'merge_results', 'output')
        os.mkdir(output_path)

        merge_files(output_path, input_path)
        set1 = self._read_file(path.join(output_path, 'subdir1', 'subsubdir1', 'file1.txt'))
        assert_that(set1, equal_to({'test1.1', 'test1.2', 'test9.1', 'test9.2'}))
        set2 = self._read_file(path.join(output_path, 'subdir1', 'subsubdir1', 'file2.txt'))
        assert_that(set2, equal_to({'test2.1', 'test2.2', 'test10.1', 'test10.2'}))
        set3 = self._read_file(path.join(output_path, 'subdir1', 'subsubdir2', 'file1.txt'))
        assert_that(set3, equal_to({'test3.1', 'test3.2', 'test11.1', 'test11.2'}))
        set4 = self._read_file(path.join(output_path, 'subdir1', 'subsubdir2', 'file2.txt'))
        assert_that(set4, equal_to({'test4.1', 'test4.2', 'test12.1', 'test12.2'}))
        set5 = self._read_file(path.join(output_path, 'subdir2', 'subsubdir1', 'file1.txt'))
        assert_that(set5, equal_to({'test5.1', 'test5.2', 'test13.1', 'test13.2'}))
        set6 = self._read_file(path.join(output_path, 'subdir2', 'subsubdir1', 'file2.txt'))
        assert_that(set6, equal_to({'test6.1', 'test6.2', 'test14.1', 'test14.2'}))
        set7 = self._read_file(path.join(output_path, 'subdir2', 'subsubdir2', 'file1.txt'))
        assert_that(set7, equal_to({'test7.1', 'test7.2', 'test15.1', 'test15.2'}))
        set8 = self._read_file(path.join(output_path, 'subdir2', 'subsubdir2', 'file2.txt'))
        assert_that(set8, equal_to({'test8.1', 'test8.2', 'test16.1', 'test16.2'}))

        shutil.rmtree(output_path)

    def _read_file(self, file_path):
        result = set()
        with open(file_path, 'r') as f:
            for line in f:
                sanitized_line = line.strip()
                if sanitized_line:
                    result.add(sanitized_line)
        return result
