import argparse
from collections import defaultdict
import numpy as np
import os
from enum import Enum
import itertools
import pandas as pd
from utils import get_all_files_paths
import matplotlib.pyplot as plt
import brewer2mpl
from scipy import stats


def agg_row(data_file_path):
    values = []
    line_no = 0
    with open(data_file_path, 'r') as data_file:
        for line in data_file:
            values.append(int(line))
            line_no += 1
            if line_no >= 100:
                break
    return values


def agg_rowcol(data_file_path):
    values = []
    line_no = 0
    with open(data_file_path, 'r') as data_file:
        for line in data_file:
            if 'NA' in line:
                print 'Found NA in file: %s' % data_file_path
                continue
            values.extend([float(v) for v in line.split(',')])
            line_no += 1
            if line_no >= 100:
                break
    #if line_no < 100:
    #    raise Exception('Found only %d proper values of 100' % line_no)
    return values


class DataType(Enum):
    found_solution_number = 'Average number of found solutions'
    eps = 'Average epsilon value'
    relations_number = 'Average number of necessary relations'
    new_relations_number = 'Average number of new necessary relations'


class MethodType(Enum):
    equal_freq = 'EF'
    equal_width = 'EW'
    ghaderi = 'GH'
    kernel = 'KDE'
    kmeans = 'KM'

    @staticmethod
    def get_ordered():
        return [
            str(MethodType.equal_freq),
            str(MethodType.equal_width),
            str(MethodType.ghaderi),
            str(MethodType.kernel),
            str(MethodType.kmeans)
        ]


class DistributionType(Enum):
    skew_normal = "skew normal"
    uniform = "uniform"

    @staticmethod
    def get_ordered():
        return [
            str(DistributionType.uniform),
            str(DistributionType.skew_normal)
        ]


class DataUnit:
    def __init__(self, data_type, distribution, crit_number, alt_number, pref_info, char_points, method_name, path):
        self.data_type = data_type
        self.distribution = distribution
        self.crit_number = crit_number
        self.alt_number = alt_number
        self.pref_info = pref_info
        self.char_points = char_points
        self.method_name = method_name
        self.path = path


class PathInterpreter:
    def __init__(self, base_path):
        self.base_path = base_path

    def interpret(self, path):
        path_parts = path.split('/')

        data_type = self.get_data_type_name(path_parts[4])
        if data_type is None:
            return None

        distribution = self.get_distribution_name(path_parts[0])
        crit_number, alt_number = [int(i) for i in path_parts[1].split('x')]
        pref_info = int(path_parts[2][len('pref-'):]) if path_parts[2].startswith('pref-')\
            else path_parts[2].lower()
        char_points = int(path_parts[3][len('SEGMENTED-')]) if path_parts[3].startswith('SEGMENTED-')\
            else path_parts[3].lower()
        method_id = path_parts[3][len('SEGMENTED-N-'):] if path_parts[3].startswith('SEGMENTED-') else None
        path = os.path.join(self.base_path, path)

        return DataUnit(data_type, distribution, crit_number, alt_number, pref_info, char_points,
                        self.get_method_name(method_id), path)

    def get_method_name(self, method_id):
        if method_id == 'EQUAL_FREQ_INTERVAL':
            return MethodType.equal_freq
        if method_id == 'EQUAL_WIDTH_INTERVAL':
            return MethodType.equal_width
        if method_id == 'GHADERI_DISCRETIZATION':
            return MethodType.ghaderi
        if method_id == 'KERNEL_DENSITY_ESTIMATION':
            return MethodType.kernel
        if method_id == 'K_MEANS':
            return MethodType.kmeans

    def get_distribution_name(self, distribution_id):
        if distribution_id == "SKEW_NORMAL":
            return DistributionType.skew_normal
        if distribution_id == "UNIFORM":
            return DistributionType.uniform

    def get_data_type_name(self, data_type_id):
        if data_type_id == 'foundsolutionsnumber.csv':
            return DataType.found_solution_number
        elif data_type_id == 'epsvalues.csv':
            return DataType.eps
        elif data_type_id == 'relationsnumbers.csv':
            return DataType.relations_number
        elif data_type_id == 'prefindrelationsnumbers.csv':
            return DataType.new_relations_number


class Aggregator:
    def __init__(self):
        self.by_crit_number = defaultdict(lambda: defaultdict(list))
        self.by_alt_number = defaultdict(lambda: defaultdict(list))
        self.by_comparison_number = defaultdict(lambda: defaultdict(list))
        self.by_distribution = defaultdict(lambda: defaultdict(list))

    def flat_data(self):
        data_list = [
            self.by_crit_number,
            self.by_alt_number,
            self.by_comparison_number,
            self.by_distribution
        ]
        for data in data_list:
            for x, x_dict in data.iteritems():
                for y, y_list in x_dict.iteritems():
                    data[x][y] = np.mean(y_list)


class DataAggregation:
    def __init__(self, output_path):
        self.output_path = output_path

    def add_data(self, data):
        raise NotImplementedError()


class AggregationByCharPoints(DataAggregation):
    def __init__(self, output_path):
        DataAggregation.__init__(self, output_path)
        self.aggregators = {
            DataType.found_solution_number: Aggregator(),
            DataType.eps: Aggregator(),
            DataType.relations_number: Aggregator(),
            DataType.new_relations_number: Aggregator()
        }

    def add_data(self, data):
        if data.method_name is None or data.method_name == MethodType.equal_width:
            aggregator = self.aggregators[data.data_type]
            agg_func = agg_row if data.data_type == DataType.found_solution_number else agg_rowcol

            ch_point_str = data.char_points if isinstance(data.char_points, str) else '%d char. p.' % data.char_points

            aggregator.by_crit_number[ch_point_str][data.crit_number].extend(agg_func(data.path))
            aggregator.by_alt_number[ch_point_str][data.alt_number].extend(agg_func(data.path))
            if data.pref_info != 'ranking':
                aggregator.by_comparison_number[ch_point_str][data.pref_info].extend(agg_func(data.path))
            aggregator.by_distribution[ch_point_str][data.distribution].extend(agg_func(data.path))

    def generate_charts(self):
        for data_type, aggregator in self.aggregators.iteritems():
            aggregator.flat_data()
            series_order = ['linear', '3 char. p.', '4 char. p.', '5 char. p.', '6 char. p.', 'general']
            self._plot_chart(aggregator.by_crit_number, x_label='Criteria number', y_label=data_type,
                             series_order=series_order)
            self._plot_chart(aggregator.by_alt_number, x_label='Alternatives number', y_label=data_type,
                             series_order=series_order)
            self._plot_chart(aggregator.by_comparison_number, x_label='Number of pairwise alternatives comparisons',
                             y_label=data_type, series_order=series_order)
            self._plot_bar_chart(aggregator.by_distribution, x_label='Distribution', y_label=data_type,
                                 series_order=series_order)

    def _plot_bar_chart(self, data, x_label, y_label, series_order=None):
        df = pd.DataFrame(data)
        if series_order:
            df = df[series_order]
        bmap = brewer2mpl.get_map('Set2', 'qualitative', len(data), reverse=True)
        colors = bmap.mpl_colors
        ax = df.plot(kind='bar', color=colors)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        locs, labels = plt.xticks()
        plt.setp(labels, rotation=0)
        plt.savefig(os.path.join(self.output_path, 'charac_point_%s_%s.pdf' % (x_label, y_label)))
        plt.close()

    def _plot_chart(self, data, x_label, y_label, set_xticks=True, series_order=None):
        df = pd.DataFrame(data)
        if series_order:
            df = df[series_order]
        bmap = brewer2mpl.get_map('Set2', 'qualitative', len(data), reverse=True)
        colors = bmap.mpl_colors
        ax = df.plot(marker='o', color=colors)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        if set_xticks:
            x_values = set()
            for series_name, series_data in data.iteritems():
                x_values.update(series_data.keys())
            ax.set_xticks(list(x_values))
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.savefig(os.path.join(self.output_path, 'charac_point_%s_%s.pdf' % (x_label, y_label)))
        plt.close()


class AggregationForMethod(DataAggregation):
    def __init__(self, output_path, method_name):
        DataAggregation.__init__(self, output_path)
        self.method_name = method_name
        self.aggregators = {
            DataType.found_solution_number: Aggregator(),
            DataType.eps: Aggregator(),
            DataType.relations_number: Aggregator(),
            DataType.new_relations_number: Aggregator()
        }

    def add_data(self, data):
        if not isinstance(data.char_points, str) and data.method_name == self.method_name:
            aggregator = self.aggregators[data.data_type]

            ch_point_str = '%d char. p.' % data.char_points
            agg_func = agg_row if data.data_type == DataType.found_solution_number else agg_rowcol

            aggregator.by_crit_number[ch_point_str][data.crit_number].extend(agg_func(data.path))
            aggregator.by_alt_number[ch_point_str][data.alt_number].extend(agg_func(data.path))
            data.pref_info = data.pref_info if data.pref_info != 'ranking' else 10
            aggregator.by_comparison_number[ch_point_str][data.pref_info].extend(agg_func(data.path))
            aggregator.by_distribution[ch_point_str][data.distribution].extend(agg_func(data.path))

    def generate_charts(self):
        for data_type, aggregator in self.aggregators.iteritems():
            aggregator.flat_data()
            xticks = [2, 4, 6, 8, 10]
            xtickslabels = ['2', '4', '6', '8', 'ranking']
            if data_type == DataType.new_relations_number or data_type == DataType.relations_number:
                xticks = [2, 4, 6, 8]
                xtickslabels = ['2', '4', '6', '8']

            self._plot_chart(aggregator.by_crit_number, x_label='Criteria number', y_label=data_type)
            self._plot_chart(aggregator.by_alt_number, x_label='Alternatives number', y_label=data_type)
            self._plot_chart(aggregator.by_comparison_number, x_label='Number of pairwise alternatives comparisons',
                             y_label=data_type, set_xticks=False, xticks=xticks, xtickslabels=xtickslabels)
            self._plot_bar_chart(aggregator.by_distribution, x_label='Distribution', y_label=data_type)

    def _plot_bar_chart(self, data, x_label, y_label):
        df = pd.DataFrame(data)
        bmap = brewer2mpl.get_map('Set2', 'qualitative', len(data), reverse=True)
        colors = bmap.mpl_colors
        ax = df.plot(kind='bar', color=colors)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        locs, labels = plt.xticks()
        plt.setp(labels, rotation=0)
        plt.savefig(os.path.join(self.output_path, '%s_%s_%s.pdf' % (self.method_name, x_label, y_label)))
        plt.close()

    def _plot_chart(self, data, x_label, y_label, set_xticks=True, xticks=None, xtickslabels=None):
        df = pd.DataFrame(data)
        bmap = brewer2mpl.get_map('Set2', 'qualitative', len(data), reverse=True)
        colors = bmap.mpl_colors
        ax = df.plot(marker='o', color=colors)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        if set_xticks:
            x_values = set()
            for series_name, series_data in data.iteritems():
                x_values.update(series_data.keys())
            ax.set_xticks(list(x_values))
        if xticks and xtickslabels:
            ax.set_xticks(xticks)
            ax.set_xticklabels(xtickslabels)
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.savefig(os.path.join(self.output_path, '%s_%s_%s.pdf' % (self.method_name, x_label, y_label)))
        plt.close()


class SummaryAggregationByMethod(DataAggregation):
    def __init__(self, output_path):
        DataAggregation.__init__(self, output_path)
        self.aggregators = {
            DataType.found_solution_number: defaultdict(lambda: defaultdict(list)),
            DataType.eps: defaultdict(lambda: defaultdict(list)),
            DataType.relations_number: defaultdict(lambda: defaultdict(list)),
            DataType.new_relations_number: defaultdict(lambda: defaultdict(list))
        }

    def add_data(self, data):
        if not isinstance(data.char_points, str):
            aggregator = self.aggregators[data.data_type]

            ch_point_str = '%d char. p.' % data.char_points
            agg_func = agg_row if data.data_type == DataType.found_solution_number else agg_rowcol

            aggregator[ch_point_str][data.method_name].extend(agg_func(data.path))

    def generate_charts(self):
        for data_type, aggregator in self.aggregators.iteritems():
            for x, x_dict in aggregator.iteritems():
                for y, y_list in x_dict.iteritems():
                    aggregator[x][y] = np.mean(y_list)
            ylim = None
            if data_type == DataType.relations_number:
                ylim = (0, 19)
            self._plot_chart(aggregator, x_label='Dicretization method code', y_label=data_type,  set_xticks=False,
                             ylim=ylim)

    def _plot_chart(self, data, x_label, y_label, set_xticks=True, xticks=None, xtickslabels=None, ylim=None):
        df = pd.DataFrame(data)
        bmap = brewer2mpl.get_map('Set2', 'qualitative', len(data), reverse=True)
        colors = bmap.mpl_colors
        ax = df.plot(kind='bar', color=colors)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        if set_xticks:
            x_values = set()
            for series_name, series_data in data.iteritems():
                x_values.update(series_data.keys())
            ax.set_xticks(list(x_values))
        if xticks and xtickslabels:
            ax.set_xticks(xticks)
            ax.set_xticklabels(xtickslabels)
        if ylim is not None:
            ax.set_ylim(ylim)
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        locs, labels = plt.xticks()
        plt.setp(labels, rotation=0)
        plt.savefig(os.path.join(self.output_path, 'summary_methods_%s_%s.pdf' % (x_label, y_label)))
        plt.close()


class WilcoxonForMethods(DataAggregation):
    class IdShortener:
        def __init__(self):
            self.artificial_ids = {}

        def shorten_id(self, identifier):
            if identifier not in self.artificial_ids:
                self.artificial_ids[identifier] = len(self.artificial_ids)
            return self.artificial_ids[identifier]

    class WilcoxonAggregator:
        def __init__(self):
            self.dist_id_shortener = WilcoxonForMethods.IdShortener()
            self.pref_id_shortener = WilcoxonForMethods.IdShortener()
            self.ch_p_id_shortener = WilcoxonForMethods.IdShortener()
            self._general_method_values = defaultdict(lambda: defaultdict(list))
            self.method_values_by_crits_no = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
            self.method_values_by_alts_no = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
            self.method_values_by_comps_no = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
            self.method_values_by_distr = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
            self.method_values_by_ch_p = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

        def agg_data(self, data, values):
            if data.method_name:
                params_code = self._get_params_code(data)
                self._general_method_values[params_code][data.method_name].extend(values)
                self.method_values_by_crits_no[data.crit_number][params_code][data.method_name].extend(values)
                self.method_values_by_alts_no[data.alt_number][params_code][data.method_name].extend(values)
                self.method_values_by_comps_no[data.pref_info][params_code][data.method_name].extend(values)
                self.method_values_by_distr[data.distribution][params_code][data.method_name].extend(values)
                self.method_values_by_ch_p[data.char_points][params_code][data.method_name].extend(values)

        def get_general_p_value(self, method1, method2):
            method1_values = []
            method2_values = []
            for values_by_method in self._general_method_values.itervalues():
                method1_values.extend(values_by_method[method1])
                method2_values.extend(values_by_method[method2])
            return stats.wilcoxon(method1_values, method2_values)[1]

        def get_p_value_for_param(self, method1, method2, values_for_param):
            p_values_by_param = {}
            for param, values_for_param in values_for_param.iteritems():
                method1_values_for_param = []
                method2_values_for_param = []
                for values_by_method in values_for_param.itervalues():
                    method1_values_for_param.extend(values_by_method[method1])
                    method2_values_for_param.extend(values_by_method[method2])
                p_value_for_param = stats.wilcoxon(method1_values_for_param, method2_values_for_param)[1]
                p_values_by_param[param] = p_value_for_param
            return p_values_by_param

        def _get_params_code(self, data):
            """
            Get unique shorten code for given parameters
            :param data: parameters
            :type data: DataUnit
            :return: int
            """
            return '%d;%d;%d;%d;%d' % (
                self.dist_id_shortener.shorten_id(data.distribution),
                data.crit_number,
                data.alt_number,
                self.pref_id_shortener.shorten_id(data.pref_info),
                self.ch_p_id_shortener.shorten_id(data.char_points)
            )

    def __init__(self, output_path):
        DataAggregation.__init__(self, output_path)
        self.aggregators = {
            DataType.found_solution_number: WilcoxonForMethods.WilcoxonAggregator(),
            DataType.eps: WilcoxonForMethods.WilcoxonAggregator(),
            DataType.relations_number: WilcoxonForMethods.WilcoxonAggregator(),
            DataType.new_relations_number: WilcoxonForMethods.WilcoxonAggregator()
        }

    def add_data(self, data):
        aggregator = self.aggregators[data.data_type]
        agg_func = agg_row if data.data_type == DataType.found_solution_number else agg_rowcol
        aggregator.agg_data(data, agg_func(data.path))

    def generate_wilcoxon_comparisons(self):
        self._save_wilcoxon_results('found_solution_number', self.aggregators[DataType.found_solution_number])
        self._save_wilcoxon_results('eps', self.aggregators[DataType.eps])
        self._save_wilcoxon_results('relations_number', self.aggregators[DataType.relations_number],
                                    robustness_exp=True)
        self._save_wilcoxon_results('new_relations_number', self.aggregators[DataType.new_relations_number],
                                    robustness_exp=True)

    def _save_wilcoxon_results(self, name, aggregator, robustness_exp=False):
        self._save_general_wilcoxon_matrix(name, aggregator)
        self._save_wilcoxon_matrix_by_param(
            '%s_criteria' % name,
            lambda m1, m2: aggregator.get_p_value_for_param(m1, m2, aggregator.method_values_by_crits_no),
            sorted(aggregator.method_values_by_crits_no))
        self._save_wilcoxon_matrix_by_param(
            '%s_alternatives' % name,
            lambda m1, m2: aggregator.get_p_value_for_param(m1, m2, aggregator.method_values_by_alts_no),
            sorted(aggregator.method_values_by_alts_no))
        comp_values = [2, 4, 6, 8, 'ranking'] if not robustness_exp else [2, 4, 6, 8]
        self._save_wilcoxon_matrix_by_param(
            '%s_comparisons_no' % name,
            lambda m1, m2: aggregator.get_p_value_for_param(m1, m2, aggregator.method_values_by_comps_no),
            comp_values)
        self._save_wilcoxon_matrix_by_param(
            '%s_distribution' % name,
            lambda m1, m2: aggregator.get_p_value_for_param(m1, m2, aggregator.method_values_by_distr),
            DistributionType.get_ordered())
        self._save_wilcoxon_matrix_by_param(
            '%s_characteristic_points' % name,
            lambda m1, m2: aggregator.get_p_value_for_param(m1, m2, aggregator.method_values_by_ch_p),
            sorted(aggregator.method_values_by_ch_p))

    def _save_general_wilcoxon_matrix(self, name, aggregator):
        with open(os.path.join(self.output_path, 'wilcoxon_%s.csv' % name), 'w') as output_file:
            methods = MethodType.get_ordered()
            output_file.write(',%s\n' % ','.join(methods))
            for method1 in methods:
                output_file.write('%s' % method1)
                for method2 in methods:
                    if method1 != method2:
                        output_file.write(',%f' % aggregator.get_general_p_value(method1, method2))
                    else:
                        output_file.write(',')
                output_file.write('\n')

    def _save_wilcoxon_matrix_by_param(self, name, p_value_for_methods_pair, parameter_values):
        with open(os.path.join(self.output_path, 'wilcoxon_%s.csv' % name), 'w') as output_file:
            methods = MethodType.get_ordered()
            header_fields = ['%s & %s' % (m1, m2) for m1, m2 in itertools.combinations(methods, 2)]
            output_file.write(',%s\n' % ','.join(header_fields))
            for param in parameter_values:
                output_file.write('%s' % str(param))
                for method1, method2 in itertools.combinations(methods, 2):
                    p_values_for_params = p_value_for_methods_pair(method1, method2)
                    if param not in p_values_for_params:
                        raise Exception("Parameter \"%s\" not found in wilcoxon results" % param)
                    output_file.write(',%f' % p_values_for_params[param])
                output_file.write('\n')


def collect_data(aggregation_method, data_path):
    files_paths = get_all_files_paths(data_path)
    interpreter = PathInterpreter(data_path)
    for file_idx, file_path in enumerate(files_paths):
        data_unit = interpreter.interpret(file_path)
        if data_unit is not None:
            aggregation_method.add_data(data_unit)
        if (file_idx + 1) % 1000 == 0:
            print 'Number of preprocessed files: %d/%d' % (file_idx + 1, len(files_paths))


def aggregate_data(output_path, data_path):
    # char_points_aggregation = AggregationByCharPoints(output_path)
    # collect_data(char_points_aggregation, data_path)
    # char_points_aggregation.generate_charts()
    #
    # method_aggregation = AggregationForMethod(output_path, MethodType.equal_width)
    # collect_data(method_aggregation, data_path)
    # method_aggregation.generate_charts()
    #
    # method_aggregation = AggregationForMethod(output_path, MethodType.equal_freq)
    # collect_data(method_aggregation, data_path)
    # method_aggregation.generate_charts()
    #
    # method_aggregation = AggregationForMethod(output_path, MethodType.ghaderi)
    # collect_data(method_aggregation, data_path)
    # method_aggregation.generate_charts()
    #
    # method_aggregation = AggregationForMethod(output_path, MethodType.kernel)
    # collect_data(method_aggregation, data_path)
    # method_aggregation.generate_charts()
    #
    # method_aggregation = AggregationForMethod(output_path, MethodType.kmeans)
    # collect_data(method_aggregation, data_path)
    # method_aggregation.generate_charts()
    #
    # summary_aggregation = SummaryAggregationByMethod(output_path)
    # collect_data(summary_aggregation, data_path)
    # summary_aggregation.generate_charts()

    wilcoxon_for_methods = WilcoxonForMethods(output_path)
    collect_data(wilcoxon_for_methods, data_path)
    wilcoxon_for_methods.generate_wilcoxon_comparisons()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("output_path", help="output path")
    parser.add_argument("data_path", help="path to merged results of experiments")
    args = parser.parse_args()

    aggregate_data(args.output_path, args.data_path)
