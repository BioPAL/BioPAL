# Python numerical libraries
import numpy as np

# Python XMl libraries
from lxml import etree as xml_tools
import xml.etree.ElementTree as ET

# Python IO
import copy

class XmlIO:
    """
    Class to read/write XML documents
    Allows python object-like access to parameters.
    pp = XmlIO('path_to_pp.xml')
    print(pp.v0)     # v0 is automatically converted to floating point
    print(pp.r[:10]) # first 10 values of range vector, also floating point
    The parameter object also allows dictionary-like access to handle problematic parameter names
    (which clash with python keywords). For example:
    print(pp['lambda']) # pp.lambda would be a syntax error
    print(pp['pass'])   # same as above
    print(pp['v0'])     # dictionary style access works for other parameters, too!
    The class provides full read/write support. Parameter values are changed by standard assignment
    and structures can be saved using the write method:
    pp.v0 = 100
    pp.write('path_to_new_pp.xml')
    """

    def __init__(self, root):
        if isinstance(root, str):
            self.__dict__['__root__'] = xml_tools.parse(root).find('object')
        else:
            self.__dict__['__root__'] = root

        if self.__root__ is None:
            raise ValueError('Expected an "object" element below the root element!')

        self.__dict__['__iterend__'] = False

    def __getstate__(self):
        return self.__root__

    def __setstate__(self, root):
        self.__dict__['__root__'] = root
        self.__dict__['__iterend__'] = False

    def __getparam__(self, name):
        p = [p for p in self.__root__.iter('parameter') if p.attrib['name'] == name]
        if len(p) != 1:
            raise AttributeError('Expected a unique match parameter name "%s", got %i matches.' % (name, len(p)))

        return [p[0].find(tag) for tag in ('remark', 'datatype', 'value', 'unit')]

    @staticmethod
    def xml2val(v, t):
        type = t.text
        shape = t.attrib['length']
        shape = np.asarray([np.uint64(d) for d in shape.split()])[::-1]
        size = np.prod(shape)

        if (type == 'pointer'):
            p = v.find('parameter')
            return XmlIO.xml2val(*[p.find(t) for t in ('value', 'datatype')])

        if (type == 'struct'):
            obj_arr = [XmlIO(obj) for obj in v.iter('object')]
            return obj_arr[0] if size <= 1 else obj_arr

        conv = {'bool': bool, 'int': int, 'long': int, 'float': np.float, 'double': np.double, 'string': lambda s: s}
        try:
            if size > 1:
                val = np.asarray([conv[type](v) for v in v.text.strip('[]').split(',')]).reshape(shape)
            else:
                val = conv[type](v.text)
        except KeyError:
            print('XmlIO WARNING: Unsupported data type "%s" encountered. Skipping!' % (type))
            return None

        return val

    @staticmethod
    def val2xml(v, t, value):
        cdict = {str: (str, 'string'),
                 int: (str, 'long'),
                 float: (str, 'double'),
                 complex: (lambda z: '({},{})'.format(z.real, z.imag), 'complex')}
        cdict[np.uint8] = cdict[int]
        cdict[np.int32] = cdict[int]
        cdict[np.uint32] = cdict[int]
        cdict[np.int64] = cdict[int]
        cdict[np.uint64] = cdict[int]
        cdict[np.float32] = (str, 'float')
        cdict[np.float64] = cdict[float]
        cdict[np.complex64] = cdict[complex]
        cdict[bool] = cdict[str]

        if (t.text == 'pointer'):
            p = v.find('parameter')
            return XmlIO.val2xml(*([p.find(t) for t in ('value', 'datatype')] + [value]))

        try:
            vsize = 1 if isinstance(value, str) else len(value)
        except TypeError:
            vsize = 1

        t.attrib['length'] = str(vsize)
        v.clear()
        if vsize == 1 and not isinstance(value, XmlIO):
            t.text = cdict[type(value)][1]
            v.text = cdict[type(value)][0](value)
        elif all([isinstance(v, XmlIO) for v in value]):
            t.text = 'struct'
            for obj in value:
                v.append(copy.deepcopy(obj.__root__))
        else:
            if isinstance(value, np.ndarray):
                t.attrib['length'] = ' '.join([str(l) for l in value.shape[::-1]])
                value = value.flat
            vtype = type(value[0])
            t.text = cdict[vtype][1]
            v.text = '[' + ', '.join([cdict[vtype][0](val) for val in value]) + ']'

    def __getattr__(self, key):
        if key in self.__dict__:
            return self.__dict__[key]
        if key == 0:
            return self
        r, t, v, u = self.__getparam__(key)
        return XmlIO.xml2val(v, t)

    def __getitem__(self, key):
        return self.__getattr__(key)

    def __setattr__(self, name, value):
        if name in self.__dict__:
            self.__dict__[name] = value
            return
        r, t, v, u = self.__getparam__(name)
        XmlIO.val2xml(v, t, value)

    def __setitem__(self, key, value):
        self.__setattr__(key, value)

    def __contains__(self, key):
        try:
            _ = self.__getparam__(key)
        except AttributeError:
            return False
        return True

    def __len__(self):
        return 1

    def __iter__(self):
        self.__iterend__ = False
        return self

    def __next__(self):
        if self.__iterend__:
            raise StopIteration()
        self.__iterend__ = True
        return self

    def update(self, obj):
        try:
            d = obj.__dict__
        except AttributeError:
            d = obj
        if not isinstance(d, dict):
            raise ValueError('Expected a dictionary or an object with a __dict__ attribute!')

        for k in d:
            self.__setattr__(k, d[k])

    def __totree(self):
        ste_root = xml_tools.Element('stexml')
        ste_root.text = '\n'
        ste_root.append(copy.deepcopy(self.__root__))
        ste_root.addprevious(xml_tools.PI('xml-stylesheet', 'type="text/xsl" href="stexml.xsl"'))
        tree = xml_tools.ElementTree(ste_root)
        return tree

    def write(self, filename):
        self.__totree().write(filename, pretty_print=True, encoding='UTF-8', xml_declaration=True)

    def tostring(self):
        return xml_tools.tostring(self.__totree().getroot(), encoding='UTF-8')

def add_param(root, name, unit_text, datatype_text, remark_text="none", value_text=None, length=1, sub_flag=False):

    # Sub flag is needed in order to correctly add parameters to a struct-type xml tree
    if sub_flag:
        # Check if the value and object subelements were created previously
        if root.find('value') is None:
            value = ET.SubElement(root, "value")
            object = ET.SubElement(value, "object")
        else:
            root = root.find('value')
            object = root.find('object')
        root = object

    param = ET.SubElement(root, "parameter")
    param.set("name", name)

    unit = ET.SubElement(param, "unit")
    unit.text = unit_text

    datatype = ET.SubElement(param, "datatype")
    datatype.set("length", str(length))
    datatype.text = datatype_text

    remark = ET.SubElement(param, "remark")
    remark.text = remark_text

    if datatype_text != 'struct':
        value = ET.SubElement(param, "value")
        value.text = value_text

    return param