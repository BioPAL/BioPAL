"""
Generic polynomial module
-------------------------
"""


class GenericPoly:
    """
    Generic polynomial class
    """
    def __init__(self, reference_values, coefficients, powers):
        """
        Initialize the GenericPoly class

        :param reference_values:
        :param coefficients:
        :param powers:
        """

        self.reference_values = reference_values
        self.poly = dict()
        for index_coefficients, coefficient in enumerate(coefficients):
            self.poly[coefficient] = powers[index_coefficients]

    def __repr__(self):
        representation = "coefficients: {}\n".format(list(self.poly.keys()))
        representation += "powers: {}\n".format(list(self.poly.values()))
        representation += "reference values: {}\n".format(list(self.reference_values))
        return representation

    def evaluate(self, values):
        """
        Evaluate the GenericPoly for the values provided

        :param values:
        :return:
        """
        values = tuple(values)
        result = 0
        for coefficient, powers in self.poly.items():
            current_result = 1
            for index_dimensions, ref_val in enumerate(self.reference_values):
                current_result *= \
                    (values[index_dimensions] - ref_val) ** powers[index_dimensions]
            result += coefficient * current_result

        return result


class SortedPolyList:
    """
    SortedPolyList class
    """
    def __init__(self, reference_index=0, list_generic_poly=None):
        """
        Initialize sorted poly list

        :param reference_index:
        :param list_generic_poly:
        """
        self.reference_index = reference_index
        if list_generic_poly is not None:
            self._sorted_poly_list = list_generic_poly
            self._sort_poly_list()
        else:
            self._sorted_poly_list = list()

    def __repr__(self):
        representation = "Sorted Poly List:\n"
        for index, poly in enumerate(self._sorted_poly_list):
            representation += "Poly #{}: \n{}\n".format(index, poly)
        return representation

    def append(self, generic_poly):
        """
        Append a new poly to the SortedPolyList

        :param generic_poly:
        """
        self._sorted_poly_list.append(generic_poly)
        self._sort_poly_list()

    def _sort_poly_list(self):
        self._sorted_poly_list.sort(key=lambda x: x.reference_values[self.reference_index])

    def evaluate(self, values):
        """
        Evaluate the SortedPolyList for the values provided

        :param values:
        :return:
        """
        previous_poly = self._sorted_poly_list[0]
        for poly in self._sorted_poly_list:
            if poly.reference_values[self.reference_index] > values[self.reference_index]:
                break
            previous_poly = poly
        return previous_poly.evaluate(values)
