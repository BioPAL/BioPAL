# arepytools/io/parsing/pyxb_metadata_types.py
# -*- coding: utf-8 -*-
# PyXB bindings for NM:09b9bbfd44fc1c40c2a0e027bb50aade954f7a67
# Generated 2021-02-17 11:26:26.605491 by PyXB version 1.2.6 using Python 3.8.5.final.0
# Namespace aresysTypes [xmlns:ns1]

from __future__ import unicode_literals
import pyxb
import pyxb.binding
import pyxb.binding.saxer
import io
import pyxb.utils.utility
import pyxb.utils.domutils
import sys
import pyxb.utils.six as _six
# Unique identifier for bindings created at the same time
_GenerationUID = pyxb.utils.utility.UniqueIdentifier('urn:uuid:972c4fe8-710a-11eb-bb6b-7085c2d17c63')

# Version of PyXB used to generate the bindings
_PyXBVersion = '1.2.6'
# Generated bindings are not compatible across PyXB versions
if pyxb.__version__ != _PyXBVersion:
    raise pyxb.PyXBVersionError(_PyXBVersion)

# A holder for module-level binding classes so we can access them from
# inside class definitions where property names may conflict.
_module_typeBindings = pyxb.utils.utility.Object()

# Import bindings for namespaces imported into schema
import pyxb.binding.datatypes

# NOTE: All namespace declarations are reserved within the binding
Namespace = pyxb.namespace.NamespaceForURI('aresysTypes', create_if_missing=True)
Namespace.configureCategories(['typeBinding', 'elementBinding'])

def CreateFromDocument (xml_text, default_namespace=None, location_base=None):
    """Parse the given XML and use the document element to create a
    Python instance.

    @param xml_text An XML document.  This should be data (Python 2
    str or Python 3 bytes), or a text (Python 2 unicode or Python 3
    str) in the L{pyxb._InputEncoding} encoding.

    @keyword default_namespace The L{pyxb.Namespace} instance to use as the
    default namespace where there is no default namespace in scope.
    If unspecified or C{None}, the namespace of the module containing
    this function will be used.

    @keyword location_base: An object to be recorded as the base of all
    L{pyxb.utils.utility.Location} instances associated with events and
    objects handled by the parser.  You might pass the URI from which
    the document was obtained.
    """

    if pyxb.XMLStyle_saxer != pyxb._XMLStyle:
        dom = pyxb.utils.domutils.StringToDOM(xml_text)
        return CreateFromDOM(dom.documentElement, default_namespace=default_namespace)
    if default_namespace is None:
        default_namespace = Namespace.fallbackNamespace()
    saxer = pyxb.binding.saxer.make_parser(fallback_namespace=default_namespace, location_base=location_base)
    handler = saxer.getContentHandler()
    xmld = xml_text
    if isinstance(xmld, _six.text_type):
        xmld = xmld.encode(pyxb._InputEncoding)
    saxer.parse(io.BytesIO(xmld))
    instance = handler.rootObject()
    return instance

def CreateFromDOM (node, default_namespace=None):
    """Create a Python instance from the given DOM node.
    The node tag must correspond to an element declaration in this module.

    @deprecated: Forcing use of DOM interface is unnecessary; use L{CreateFromDocument}."""
    if default_namespace is None:
        default_namespace = Namespace.fallbackNamespace()
    return pyxb.binding.basis.element.AnyCreateFromDOM(node, default_namespace)


# Atomic simple type: {aresysTypes}CellTypeVerboseType
class CellTypeVerboseType (pyxb.binding.datatypes.string, pyxb.binding.basis.enumeration_mixin):

    """An atomic simple type."""

    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'CellTypeVerboseType')
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 3, 2)
    _Documentation = None
CellTypeVerboseType._CF_enumeration = pyxb.binding.facets.CF_enumeration(value_datatype=CellTypeVerboseType, enum_prefix=None)
CellTypeVerboseType.FLOAT_COMPLEX = CellTypeVerboseType._CF_enumeration.addEnumeration(unicode_value='FLOAT_COMPLEX', tag='FLOAT_COMPLEX')
CellTypeVerboseType.FLOAT32 = CellTypeVerboseType._CF_enumeration.addEnumeration(unicode_value='FLOAT32', tag='FLOAT32')
CellTypeVerboseType.DOUBLE_COMPLEX = CellTypeVerboseType._CF_enumeration.addEnumeration(unicode_value='DOUBLE_COMPLEX', tag='DOUBLE_COMPLEX')
CellTypeVerboseType.FLOAT64 = CellTypeVerboseType._CF_enumeration.addEnumeration(unicode_value='FLOAT64', tag='FLOAT64')
CellTypeVerboseType.INT16 = CellTypeVerboseType._CF_enumeration.addEnumeration(unicode_value='INT16', tag='INT16')
CellTypeVerboseType.SHORT_COMPLEX = CellTypeVerboseType._CF_enumeration.addEnumeration(unicode_value='SHORT_COMPLEX', tag='SHORT_COMPLEX')
CellTypeVerboseType.INT32 = CellTypeVerboseType._CF_enumeration.addEnumeration(unicode_value='INT32', tag='INT32')
CellTypeVerboseType.INT_COMPLEX = CellTypeVerboseType._CF_enumeration.addEnumeration(unicode_value='INT_COMPLEX', tag='INT_COMPLEX')
CellTypeVerboseType.INT8 = CellTypeVerboseType._CF_enumeration.addEnumeration(unicode_value='INT8', tag='INT8')
CellTypeVerboseType.INT8_COMPLEX = CellTypeVerboseType._CF_enumeration.addEnumeration(unicode_value='INT8_COMPLEX', tag='INT8_COMPLEX')
CellTypeVerboseType.CUSTOM = CellTypeVerboseType._CF_enumeration.addEnumeration(unicode_value='CUSTOM', tag='CUSTOM')
CellTypeVerboseType._InitializeFacetMap(CellTypeVerboseType._CF_enumeration)
Namespace.addCategoryObject('typeBinding', 'CellTypeVerboseType', CellTypeVerboseType)
_module_typeBindings.CellTypeVerboseType = CellTypeVerboseType

# Atomic simple type: {aresysTypes}RasterFormatType
class RasterFormatType (pyxb.binding.datatypes.string, pyxb.binding.basis.enumeration_mixin):

    """An atomic simple type."""

    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'RasterFormatType')
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 18, 2)
    _Documentation = None
RasterFormatType._CF_enumeration = pyxb.binding.facets.CF_enumeration(value_datatype=RasterFormatType, enum_prefix=None)
RasterFormatType.ARESYS_RASTER = RasterFormatType._CF_enumeration.addEnumeration(unicode_value='ARESYS_RASTER', tag='ARESYS_RASTER')
RasterFormatType.ARESYS_GEOTIFF = RasterFormatType._CF_enumeration.addEnumeration(unicode_value='ARESYS_GEOTIFF', tag='ARESYS_GEOTIFF')
RasterFormatType.RASTER = RasterFormatType._CF_enumeration.addEnumeration(unicode_value='RASTER', tag='RASTER')
RasterFormatType._InitializeFacetMap(RasterFormatType._CF_enumeration)
Namespace.addCategoryObject('typeBinding', 'RasterFormatType', RasterFormatType)
_module_typeBindings.RasterFormatType = RasterFormatType

# Atomic simple type: {aresysTypes}SensorNamesType
class SensorNamesType (pyxb.binding.datatypes.string, pyxb.binding.basis.enumeration_mixin):

    """An atomic simple type."""

    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'SensorNamesType')
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 25, 2)
    _Documentation = None
SensorNamesType._CF_enumeration = pyxb.binding.facets.CF_enumeration(value_datatype=SensorNamesType, enum_prefix=None)
SensorNamesType.NOT_SET = SensorNamesType._CF_enumeration.addEnumeration(unicode_value='NOT SET', tag='NOT_SET')
SensorNamesType.ASAR = SensorNamesType._CF_enumeration.addEnumeration(unicode_value='ASAR', tag='ASAR')
SensorNamesType.PALSAR = SensorNamesType._CF_enumeration.addEnumeration(unicode_value='PALSAR', tag='PALSAR')
SensorNamesType.ERS1 = SensorNamesType._CF_enumeration.addEnumeration(unicode_value='ERS1', tag='ERS1')
SensorNamesType.ERS2 = SensorNamesType._CF_enumeration.addEnumeration(unicode_value='ERS2', tag='ERS2')
SensorNamesType.RADARSAT = SensorNamesType._CF_enumeration.addEnumeration(unicode_value='RADARSAT', tag='RADARSAT')
SensorNamesType.TERRASARX = SensorNamesType._CF_enumeration.addEnumeration(unicode_value='TERRASARX', tag='TERRASARX')
SensorNamesType.SENTINEL1 = SensorNamesType._CF_enumeration.addEnumeration(unicode_value='SENTINEL1', tag='SENTINEL1')
SensorNamesType.SENTINEL1A = SensorNamesType._CF_enumeration.addEnumeration(unicode_value='SENTINEL1A', tag='SENTINEL1A')
SensorNamesType.SENTINEL1B = SensorNamesType._CF_enumeration.addEnumeration(unicode_value='SENTINEL1B', tag='SENTINEL1B')
SensorNamesType.SAOCOM = SensorNamesType._CF_enumeration.addEnumeration(unicode_value='SAOCOM', tag='SAOCOM')
SensorNamesType.SAOCOM_1A = SensorNamesType._CF_enumeration.addEnumeration(unicode_value='SAOCOM-1A', tag='SAOCOM_1A')
SensorNamesType.SAOCOM_1B = SensorNamesType._CF_enumeration.addEnumeration(unicode_value='SAOCOM-1B', tag='SAOCOM_1B')
SensorNamesType.UAVSAR = SensorNamesType._CF_enumeration.addEnumeration(unicode_value='UAVSAR', tag='UAVSAR')
SensorNamesType._InitializeFacetMap(SensorNamesType._CF_enumeration)
Namespace.addCategoryObject('typeBinding', 'SensorNamesType', SensorNamesType)
_module_typeBindings.SensorNamesType = SensorNamesType

# Atomic simple type: {aresysTypes}AcquisitionModeType
class AcquisitionModeType (pyxb.binding.datatypes.string, pyxb.binding.basis.enumeration_mixin):

    """An atomic simple type."""

    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'AcquisitionModeType')
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 43, 2)
    _Documentation = None
AcquisitionModeType._CF_enumeration = pyxb.binding.facets.CF_enumeration(value_datatype=AcquisitionModeType, enum_prefix=None)
AcquisitionModeType.NOT_SET = AcquisitionModeType._CF_enumeration.addEnumeration(unicode_value='NOT SET', tag='NOT_SET')
AcquisitionModeType.STRIPMAP = AcquisitionModeType._CF_enumeration.addEnumeration(unicode_value='STRIPMAP', tag='STRIPMAP')
AcquisitionModeType.DOUBLE_POL = AcquisitionModeType._CF_enumeration.addEnumeration(unicode_value='DOUBLE POL', tag='DOUBLE_POL')
AcquisitionModeType.QUAD_POL = AcquisitionModeType._CF_enumeration.addEnumeration(unicode_value='QUAD POL', tag='QUAD_POL')
AcquisitionModeType.SCANSAR = AcquisitionModeType._CF_enumeration.addEnumeration(unicode_value='SCANSAR', tag='SCANSAR')
AcquisitionModeType.TOPSAR = AcquisitionModeType._CF_enumeration.addEnumeration(unicode_value='TOPSAR', tag='TOPSAR')
AcquisitionModeType.SPOT = AcquisitionModeType._CF_enumeration.addEnumeration(unicode_value='SPOT', tag='SPOT')
AcquisitionModeType.WAVE = AcquisitionModeType._CF_enumeration.addEnumeration(unicode_value='WAVE', tag='WAVE')
AcquisitionModeType.GMTI = AcquisitionModeType._CF_enumeration.addEnumeration(unicode_value='GMTI', tag='GMTI')
AcquisitionModeType._InitializeFacetMap(AcquisitionModeType._CF_enumeration)
Namespace.addCategoryObject('typeBinding', 'AcquisitionModeType', AcquisitionModeType)
_module_typeBindings.AcquisitionModeType = AcquisitionModeType

# Atomic simple type: {aresysTypes}Endianity
class Endianity (pyxb.binding.datatypes.string, pyxb.binding.basis.enumeration_mixin):

    """An atomic simple type."""

    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'Endianity')
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 56, 2)
    _Documentation = None
Endianity._CF_enumeration = pyxb.binding.facets.CF_enumeration(value_datatype=Endianity, enum_prefix=None)
Endianity.BIGENDIAN = Endianity._CF_enumeration.addEnumeration(unicode_value='BIGENDIAN', tag='BIGENDIAN')
Endianity.LITTLEENDIAN = Endianity._CF_enumeration.addEnumeration(unicode_value='LITTLEENDIAN', tag='LITTLEENDIAN')
Endianity._InitializeFacetMap(Endianity._CF_enumeration)
Namespace.addCategoryObject('typeBinding', 'Endianity', Endianity)
_module_typeBindings.Endianity = Endianity

# Atomic simple type: {aresysTypes}AscendingDescendingType
class AscendingDescendingType (pyxb.binding.datatypes.string, pyxb.binding.basis.enumeration_mixin):

    """An atomic simple type."""

    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'AscendingDescendingType')
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 66, 2)
    _Documentation = None
AscendingDescendingType._CF_enumeration = pyxb.binding.facets.CF_enumeration(value_datatype=AscendingDescendingType, enum_prefix=None)
AscendingDescendingType.ASCENDING = AscendingDescendingType._CF_enumeration.addEnumeration(unicode_value='ASCENDING', tag='ASCENDING')
AscendingDescendingType.DESCENDING = AscendingDescendingType._CF_enumeration.addEnumeration(unicode_value='DESCENDING', tag='DESCENDING')
AscendingDescendingType.NOT_AVAILABLE = AscendingDescendingType._CF_enumeration.addEnumeration(unicode_value='NOT_AVAILABLE', tag='NOT_AVAILABLE')
AscendingDescendingType._InitializeFacetMap(AscendingDescendingType._CF_enumeration)
Namespace.addCategoryObject('typeBinding', 'AscendingDescendingType', AscendingDescendingType)
_module_typeBindings.AscendingDescendingType = AscendingDescendingType

# Atomic simple type: {aresysTypes}units
class units (pyxb.binding.datatypes.string, pyxb.binding.basis.enumeration_mixin):

    """An atomic simple type."""

    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'units')
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 73, 2)
    _Documentation = None
units._CF_enumeration = pyxb.binding.facets.CF_enumeration(value_datatype=units, enum_prefix=None)
units.emptyString = units._CF_enumeration.addEnumeration(unicode_value='', tag='emptyString')
units.m = units._CF_enumeration.addEnumeration(unicode_value='m', tag='m')
units.s = units._CF_enumeration.addEnumeration(unicode_value='s', tag='s')
units.j = units._CF_enumeration.addEnumeration(unicode_value='j', tag='j')
units.dB = units._CF_enumeration.addEnumeration(unicode_value='dB', tag='dB')
units.rad = units._CF_enumeration.addEnumeration(unicode_value='rad', tag='rad')
units.deg = units._CF_enumeration.addEnumeration(unicode_value='deg', tag='deg')
units.ms = units._CF_enumeration.addEnumeration(unicode_value='m/s', tag='ms')
units.ms2 = units._CF_enumeration.addEnumeration(unicode_value='m/s2', tag='ms2')
units.ms3 = units._CF_enumeration.addEnumeration(unicode_value='m/s3', tag='ms3')
units.ms4 = units._CF_enumeration.addEnumeration(unicode_value='m/s4', tag='ms4')
units.ss = units._CF_enumeration.addEnumeration(unicode_value='s/s', tag='ss')
units.ss2 = units._CF_enumeration.addEnumeration(unicode_value='s/s2', tag='ss2')
units.ss3 = units._CF_enumeration.addEnumeration(unicode_value='s/s3', tag='ss3')
units.ss4 = units._CF_enumeration.addEnumeration(unicode_value='s/s4', tag='ss4')
units.ss5 = units._CF_enumeration.addEnumeration(unicode_value='s/s5', tag='ss5')
units.Hzs = units._CF_enumeration.addEnumeration(unicode_value='Hz/s', tag='Hzs')
units.Hzs2 = units._CF_enumeration.addEnumeration(unicode_value='Hz/s2', tag='Hzs2')
units.Hzs3 = units._CF_enumeration.addEnumeration(unicode_value='Hz/s3', tag='Hzs3')
units.Hzs4 = units._CF_enumeration.addEnumeration(unicode_value='Hz/s4', tag='Hzs4')
units.Hzs5 = units._CF_enumeration.addEnumeration(unicode_value='Hz/s5', tag='Hzs5')
units.rads = units._CF_enumeration.addEnumeration(unicode_value='rad/s', tag='rads')
units.rads2 = units._CF_enumeration.addEnumeration(unicode_value='rad/s2', tag='rads2')
units.rads3 = units._CF_enumeration.addEnumeration(unicode_value='rad/s3', tag='rads3')
units.rads4 = units._CF_enumeration.addEnumeration(unicode_value='rad/s4', tag='rads4')
units.rads5 = units._CF_enumeration.addEnumeration(unicode_value='rad/s5', tag='rads5')
units.s85 = units._CF_enumeration.addEnumeration(unicode_value='s85', tag='s85')
units.Utc = units._CF_enumeration.addEnumeration(unicode_value='Utc', tag='Utc')
units.b = units._CF_enumeration.addEnumeration(unicode_value='b', tag='b')
units.Hz = units._CF_enumeration.addEnumeration(unicode_value='Hz', tag='Hz')
units.K = units._CF_enumeration.addEnumeration(unicode_value='K', tag='K')
units.sm = units._CF_enumeration.addEnumeration(unicode_value='s/m', tag='sm')
units.sm2 = units._CF_enumeration.addEnumeration(unicode_value='s/m2', tag='sm2')
units.sm3 = units._CF_enumeration.addEnumeration(unicode_value='s/m3', tag='sm3')
units.sm4 = units._CF_enumeration.addEnumeration(unicode_value='s/m4', tag='sm4')
units.degs = units._CF_enumeration.addEnumeration(unicode_value='deg/s', tag='degs')
units.degs2 = units._CF_enumeration.addEnumeration(unicode_value='deg/s2', tag='degs2')
units.degs3 = units._CF_enumeration.addEnumeration(unicode_value='deg/s3', tag='degs3')
units.degs4 = units._CF_enumeration.addEnumeration(unicode_value='deg/s4', tag='degs4')
units._InitializeFacetMap(units._CF_enumeration)
Namespace.addCategoryObject('typeBinding', 'units', units)
_module_typeBindings.units = units

# Atomic simple type: {aresysTypes}PolarizationType
class PolarizationType (pyxb.binding.datatypes.string, pyxb.binding.basis.enumeration_mixin):

    """An atomic simple type."""

    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'PolarizationType')
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 116, 2)
    _Documentation = None
PolarizationType._CF_enumeration = pyxb.binding.facets.CF_enumeration(value_datatype=PolarizationType, enum_prefix=None)
PolarizationType.HH = PolarizationType._CF_enumeration.addEnumeration(unicode_value='H/H', tag='HH')
PolarizationType.HV = PolarizationType._CF_enumeration.addEnumeration(unicode_value='H/V', tag='HV')
PolarizationType.VH = PolarizationType._CF_enumeration.addEnumeration(unicode_value='V/H', tag='VH')
PolarizationType.VV = PolarizationType._CF_enumeration.addEnumeration(unicode_value='V/V', tag='VV')
PolarizationType.XX = PolarizationType._CF_enumeration.addEnumeration(unicode_value='X/X', tag='XX')
PolarizationType._InitializeFacetMap(PolarizationType._CF_enumeration)
Namespace.addCategoryObject('typeBinding', 'PolarizationType', PolarizationType)
_module_typeBindings.PolarizationType = PolarizationType

# Atomic simple type: {aresysTypes}GlobalPolarizationType
class GlobalPolarizationType (pyxb.binding.datatypes.string, pyxb.binding.basis.enumeration_mixin):

    """An atomic simple type."""

    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'GlobalPolarizationType')
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 125, 2)
    _Documentation = None
GlobalPolarizationType._CF_enumeration = pyxb.binding.facets.CF_enumeration(value_datatype=GlobalPolarizationType, enum_prefix=None)
GlobalPolarizationType.SINGLE_POL = GlobalPolarizationType._CF_enumeration.addEnumeration(unicode_value='SINGLE POL', tag='SINGLE_POL')
GlobalPolarizationType.DUAL_POL = GlobalPolarizationType._CF_enumeration.addEnumeration(unicode_value='DUAL POL', tag='DUAL_POL')
GlobalPolarizationType._InitializeFacetMap(GlobalPolarizationType._CF_enumeration)
Namespace.addCategoryObject('typeBinding', 'GlobalPolarizationType', GlobalPolarizationType)
_module_typeBindings.GlobalPolarizationType = GlobalPolarizationType

# Atomic simple type: {aresysTypes}LeftRightType
class LeftRightType (pyxb.binding.datatypes.string, pyxb.binding.basis.enumeration_mixin):

    """An atomic simple type."""

    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'LeftRightType')
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 131, 2)
    _Documentation = None
LeftRightType._CF_enumeration = pyxb.binding.facets.CF_enumeration(value_datatype=LeftRightType, enum_prefix=None)
LeftRightType.LEFT = LeftRightType._CF_enumeration.addEnumeration(unicode_value='LEFT', tag='LEFT')
LeftRightType.RIGHT = LeftRightType._CF_enumeration.addEnumeration(unicode_value='RIGHT', tag='RIGHT')
LeftRightType._InitializeFacetMap(LeftRightType._CF_enumeration)
Namespace.addCategoryObject('typeBinding', 'LeftRightType', LeftRightType)
_module_typeBindings.LeftRightType = LeftRightType

# Atomic simple type: {aresysTypes}ReferenceFrameType
class ReferenceFrameType (pyxb.binding.datatypes.string, pyxb.binding.basis.enumeration_mixin):

    """An atomic simple type."""

    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'ReferenceFrameType')
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 137, 2)
    _Documentation = None
ReferenceFrameType._CF_enumeration = pyxb.binding.facets.CF_enumeration(value_datatype=ReferenceFrameType, enum_prefix=None)
ReferenceFrameType.GEOCENTRIC = ReferenceFrameType._CF_enumeration.addEnumeration(unicode_value='GEOCENTRIC', tag='GEOCENTRIC')
ReferenceFrameType.GEODETIC = ReferenceFrameType._CF_enumeration.addEnumeration(unicode_value='GEODETIC', tag='GEODETIC')
ReferenceFrameType.ZERODOPPLER = ReferenceFrameType._CF_enumeration.addEnumeration(unicode_value='ZERODOPPLER', tag='ZERODOPPLER')
ReferenceFrameType._InitializeFacetMap(ReferenceFrameType._CF_enumeration)
Namespace.addCategoryObject('typeBinding', 'ReferenceFrameType', ReferenceFrameType)
_module_typeBindings.ReferenceFrameType = ReferenceFrameType

# Atomic simple type: {aresysTypes}RotationOrderType
class RotationOrderType (pyxb.binding.datatypes.string, pyxb.binding.basis.enumeration_mixin):

    """An atomic simple type."""

    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'RotationOrderType')
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 144, 2)
    _Documentation = None
RotationOrderType._CF_enumeration = pyxb.binding.facets.CF_enumeration(value_datatype=RotationOrderType, enum_prefix=None)
RotationOrderType.YPR = RotationOrderType._CF_enumeration.addEnumeration(unicode_value='YPR', tag='YPR')
RotationOrderType.YRP = RotationOrderType._CF_enumeration.addEnumeration(unicode_value='YRP', tag='YRP')
RotationOrderType.PRY = RotationOrderType._CF_enumeration.addEnumeration(unicode_value='PRY', tag='PRY')
RotationOrderType.PYR = RotationOrderType._CF_enumeration.addEnumeration(unicode_value='PYR', tag='PYR')
RotationOrderType.RPY = RotationOrderType._CF_enumeration.addEnumeration(unicode_value='RPY', tag='RPY')
RotationOrderType.RYP = RotationOrderType._CF_enumeration.addEnumeration(unicode_value='RYP', tag='RYP')
RotationOrderType._InitializeFacetMap(RotationOrderType._CF_enumeration)
Namespace.addCategoryObject('typeBinding', 'RotationOrderType', RotationOrderType)
_module_typeBindings.RotationOrderType = RotationOrderType

# Atomic simple type: {aresysTypes}AttitudeType
class AttitudeType (pyxb.binding.datatypes.string, pyxb.binding.basis.enumeration_mixin):

    """An atomic simple type."""

    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'AttitudeType')
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 154, 2)
    _Documentation = None
AttitudeType._CF_enumeration = pyxb.binding.facets.CF_enumeration(value_datatype=AttitudeType, enum_prefix=None)
AttitudeType.NOMINAL = AttitudeType._CF_enumeration.addEnumeration(unicode_value='NOMINAL', tag='NOMINAL')
AttitudeType.REFINED = AttitudeType._CF_enumeration.addEnumeration(unicode_value='REFINED', tag='REFINED')
AttitudeType._InitializeFacetMap(AttitudeType._CF_enumeration)
Namespace.addCategoryObject('typeBinding', 'AttitudeType', AttitudeType)
_module_typeBindings.AttitudeType = AttitudeType

# Atomic simple type: [anonymous]
class STD_ANON (pyxb.binding.datatypes.string, pyxb.binding.basis.enumeration_mixin):

    """An atomic simple type."""

    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 234, 12)
    _Documentation = None
STD_ANON._CF_enumeration = pyxb.binding.facets.CF_enumeration(value_datatype=STD_ANON, enum_prefix=None)
STD_ANON.UP = STD_ANON._CF_enumeration.addEnumeration(unicode_value='UP', tag='UP')
STD_ANON.DOWN = STD_ANON._CF_enumeration.addEnumeration(unicode_value='DOWN', tag='DOWN')
STD_ANON._InitializeFacetMap(STD_ANON._CF_enumeration)
_module_typeBindings.STD_ANON = STD_ANON

# Atomic simple type: [anonymous]
class STD_ANON_ (pyxb.binding.datatypes.integer):

    """An atomic simple type."""

    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 428, 12)
    _Documentation = None
STD_ANON_._CF_minInclusive = pyxb.binding.facets.CF_minInclusive(value_datatype=STD_ANON_, value=pyxb.binding.datatypes.integer(0))
STD_ANON_._CF_maxInclusive = pyxb.binding.facets.CF_maxInclusive(value_datatype=STD_ANON_, value=pyxb.binding.datatypes.integer(99999999))
STD_ANON_._InitializeFacetMap(STD_ANON_._CF_minInclusive,
   STD_ANON_._CF_maxInclusive)
_module_typeBindings.STD_ANON_ = STD_ANON_

# Complex type {aresysTypes}TreeElementBaseType with content type EMPTY
class TreeElementBaseType (pyxb.binding.basis.complexTypeDefinition):
    """Complex type {aresysTypes}TreeElementBaseType with content type EMPTY"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_EMPTY
    _Abstract = False
    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'TreeElementBaseType')
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 62, 2)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Attribute Number uses Python identifier Number
    __Number = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'Number'), 'Number', '__aresysTypes_TreeElementBaseType_Number', pyxb.binding.datatypes.unsignedInt)
    __Number._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 63, 4)
    __Number._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 63, 4)
    
    Number = property(__Number.value, __Number.set, None, None)

    
    # Attribute Total uses Python identifier Total
    __Total = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'Total'), 'Total', '__aresysTypes_TreeElementBaseType_Total', pyxb.binding.datatypes.unsignedInt)
    __Total._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 64, 4)
    __Total._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 64, 4)
    
    Total = property(__Total.value, __Total.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __Number.name() : __Number,
        __Total.name() : __Total
    })
_module_typeBindings.TreeElementBaseType = TreeElementBaseType
Namespace.addCategoryObject('typeBinding', 'TreeElementBaseType', TreeElementBaseType)


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON (pyxb.binding.basis.complexTypeDefinition):
    """Complex type [anonymous] with content type SIMPLE"""
    _TypeDefinition = pyxb.binding.datatypes.unsignedInt
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 206, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.unsignedInt
    _ElementMap.update({
        
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.CTD_ANON = CTD_ANON


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_ (pyxb.binding.basis.complexTypeDefinition):
    """Complex type [anonymous] with content type SIMPLE"""
    _TypeDefinition = pyxb.binding.datatypes.unsignedInt
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 213, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.unsignedInt
    _ElementMap.update({
        
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.CTD_ANON_ = CTD_ANON_


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_2 (pyxb.binding.basis.complexTypeDefinition):
    """Complex type [anonymous] with content type SIMPLE"""
    _TypeDefinition = pyxb.binding.datatypes.string
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 348, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.string
    
    # Attribute Format uses Python identifier Format
    __Format = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'Format'), 'Format', '__aresysTypes_CTD_ANON_2_Format', pyxb.binding.datatypes.string, required=True)
    __Format._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 351, 18)
    __Format._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 351, 18)
    
    Format = property(__Format.value, __Format.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __Format.name() : __Format
    })
_module_typeBindings.CTD_ANON_2 = CTD_ANON_2


# Complex type [anonymous] with content type ELEMENT_ONLY
class CTD_ANON_3 (pyxb.binding.basis.complexTypeDefinition):
    """Orbit state vectors position coordinates (xyz) [m]"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 470, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Element val uses Python identifier val
    __val = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'val'), 'val', '__aresysTypes_CTD_ANON_3_val', True, pyxb.utils.utility.Location('aresysTypes.xsd', 472, 16), )

    
    val = property(__val.value, __val.set, None, None)

    _ElementMap.update({
        __val.name() : __val
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.CTD_ANON_3 = CTD_ANON_3


# Complex type [anonymous] with content type ELEMENT_ONLY
class CTD_ANON_4 (pyxb.binding.basis.complexTypeDefinition):
    """Orbit state vectors velocity coordinates [m/s]"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 489, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Element val uses Python identifier val
    __val = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'val'), 'val', '__aresysTypes_CTD_ANON_4_val', True, pyxb.utils.utility.Location('aresysTypes.xsd', 491, 16), )

    
    val = property(__val.value, __val.set, None, None)

    _ElementMap.update({
        __val.name() : __val
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.CTD_ANON_4 = CTD_ANON_4


# Complex type [anonymous] with content type ELEMENT_ONLY
class CTD_ANON_5 (pyxb.binding.basis.complexTypeDefinition):
    """Coordinates of the ascending node"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 542, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Element val uses Python identifier val
    __val = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'val'), 'val', '__aresysTypes_CTD_ANON_5_val', True, pyxb.utils.utility.Location('aresysTypes.xsd', 544, 16), )

    
    val = property(__val.value, __val.set, None, None)

    _ElementMap.update({
        __val.name() : __val
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.CTD_ANON_5 = CTD_ANON_5


# Complex type [anonymous] with content type ELEMENT_ONLY
class CTD_ANON_6 (pyxb.binding.basis.complexTypeDefinition):
    """Yaw angle values [deg]"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 601, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Element val uses Python identifier val
    __val = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'val'), 'val', '__aresysTypes_CTD_ANON_6_val', True, pyxb.utils.utility.Location('aresysTypes.xsd', 603, 16), )

    
    val = property(__val.value, __val.set, None, None)

    _ElementMap.update({
        __val.name() : __val
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.CTD_ANON_6 = CTD_ANON_6


# Complex type [anonymous] with content type ELEMENT_ONLY
class CTD_ANON_7 (pyxb.binding.basis.complexTypeDefinition):
    """Pitch angle values [deg]"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 620, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Element val uses Python identifier val
    __val = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'val'), 'val', '__aresysTypes_CTD_ANON_7_val', True, pyxb.utils.utility.Location('aresysTypes.xsd', 622, 16), )

    
    val = property(__val.value, __val.set, None, None)

    _ElementMap.update({
        __val.name() : __val
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.CTD_ANON_7 = CTD_ANON_7


# Complex type [anonymous] with content type ELEMENT_ONLY
class CTD_ANON_8 (pyxb.binding.basis.complexTypeDefinition):
    """Roll angle values [deg]"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 639, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Element val uses Python identifier val
    __val = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'val'), 'val', '__aresysTypes_CTD_ANON_8_val', True, pyxb.utils.utility.Location('aresysTypes.xsd', 641, 16), )

    
    val = property(__val.value, __val.set, None, None)

    _ElementMap.update({
        __val.name() : __val
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.CTD_ANON_8 = CTD_ANON_8


# Complex type [anonymous] with content type ELEMENT_ONLY
class CTD_ANON_9 (pyxb.binding.basis.complexTypeDefinition):
    """Azimuth antenna steering rate polynomial coefficients: const [rad/s], az [rad/s^2], az^2 [rad/s^3]"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 754, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Element val uses Python identifier val
    __val = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'val'), 'val', '__aresysTypes_CTD_ANON_9_val', True, pyxb.utils.utility.Location('aresysTypes.xsd', 756, 16), )

    
    val = property(__val.value, __val.set, None, None)

    _ElementMap.update({
        __val.name() : __val
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.CTD_ANON_9 = CTD_ANON_9


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_10 (pyxb.binding.basis.complexTypeDefinition):
    """Number of missing lines"""
    _TypeDefinition = pyxb.binding.datatypes.unsignedInt
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 804, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.unsignedInt
    _ElementMap.update({
        
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.CTD_ANON_10 = CTD_ANON_10


# Complex type [anonymous] with content type ELEMENT_ONLY
class CTD_ANON_11 (pyxb.binding.basis.complexTypeDefinition):
    """Azimuth relative times for each missing line"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 814, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Element val uses Python identifier val
    __val = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'val'), 'val', '__aresysTypes_CTD_ANON_11_val', True, pyxb.utils.utility.Location('aresysTypes.xsd', 816, 16), )

    
    val = property(__val.value, __val.set, None, None)

    _ElementMap.update({
        __val.name() : __val
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.CTD_ANON_11 = CTD_ANON_11


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_12 (pyxb.binding.basis.complexTypeDefinition):
    """Number of duplicated lines"""
    _TypeDefinition = pyxb.binding.datatypes.unsignedInt
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 833, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.unsignedInt
    _ElementMap.update({
        
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.CTD_ANON_12 = CTD_ANON_12


# Complex type [anonymous] with content type ELEMENT_ONLY
class CTD_ANON_13 (pyxb.binding.basis.complexTypeDefinition):
    """Azimuth relative times for each duplicated line"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 843, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Element val uses Python identifier val
    __val = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'val'), 'val', '__aresysTypes_CTD_ANON_13_val', True, pyxb.utils.utility.Location('aresysTypes.xsd', 845, 16), )

    
    val = property(__val.value, __val.set, None, None)

    _ElementMap.update({
        __val.name() : __val
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.CTD_ANON_13 = CTD_ANON_13


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_14 (pyxb.binding.basis.complexTypeDefinition):
    """Number of PRF changes"""
    _TypeDefinition = pyxb.binding.datatypes.unsignedInt
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 862, 14)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.unsignedInt
    _ElementMap.update({
        
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.CTD_ANON_14 = CTD_ANON_14


# Complex type [anonymous] with content type ELEMENT_ONLY
class CTD_ANON_15 (pyxb.binding.basis.complexTypeDefinition):
    """Azimuth relative times for each PRF change"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 872, 14)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Element val uses Python identifier val
    __val = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'val'), 'val', '__aresysTypes_CTD_ANON_15_val', True, pyxb.utils.utility.Location('aresysTypes.xsd', 874, 22), )

    
    val = property(__val.value, __val.set, None, None)

    _ElementMap.update({
        __val.name() : __val
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.CTD_ANON_15 = CTD_ANON_15


# Complex type [anonymous] with content type ELEMENT_ONLY
class CTD_ANON_16 (pyxb.binding.basis.complexTypeDefinition):
    """PRF changes values"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 891, 14)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Element val uses Python identifier val
    __val = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'val'), 'val', '__aresysTypes_CTD_ANON_16_val', True, pyxb.utils.utility.Location('aresysTypes.xsd', 893, 22), )

    
    val = property(__val.value, __val.set, None, None)

    _ElementMap.update({
        __val.name() : __val
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.CTD_ANON_16 = CTD_ANON_16


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_17 (pyxb.binding.basis.complexTypeDefinition):
    """Number of SWST changes"""
    _TypeDefinition = pyxb.binding.datatypes.unsignedInt
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 910, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.unsignedInt
    _ElementMap.update({
        
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.CTD_ANON_17 = CTD_ANON_17


# Complex type [anonymous] with content type ELEMENT_ONLY
class CTD_ANON_18 (pyxb.binding.basis.complexTypeDefinition):
    """Azimuth relative times for each SWST change"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 920, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Element val uses Python identifier val
    __val = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'val'), 'val', '__aresysTypes_CTD_ANON_18_val', True, pyxb.utils.utility.Location('aresysTypes.xsd', 922, 16), )

    
    val = property(__val.value, __val.set, None, None)

    _ElementMap.update({
        __val.name() : __val
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.CTD_ANON_18 = CTD_ANON_18


# Complex type [anonymous] with content type ELEMENT_ONLY
class CTD_ANON_19 (pyxb.binding.basis.complexTypeDefinition):
    """SWST changes values"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 939, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Element val uses Python identifier val
    __val = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'val'), 'val', '__aresysTypes_CTD_ANON_19_val', True, pyxb.utils.utility.Location('aresysTypes.xsd', 941, 16), )

    
    val = property(__val.value, __val.set, None, None)

    _ElementMap.update({
        __val.name() : __val
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.CTD_ANON_19 = CTD_ANON_19


# Complex type [anonymous] with content type ELEMENT_ONLY
class CTD_ANON_20 (pyxb.binding.basis.complexTypeDefinition):
    """Azimuth relative times for each noise packet"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 963, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Element val uses Python identifier val
    __val = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'val'), 'val', '__aresysTypes_CTD_ANON_20_val', True, pyxb.utils.utility.Location('aresysTypes.xsd', 965, 16), )

    
    val = property(__val.value, __val.set, None, None)

    _ElementMap.update({
        __val.name() : __val
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.CTD_ANON_20 = CTD_ANON_20


# Complex type [anonymous] with content type ELEMENT_ONLY
class CTD_ANON_21 (pyxb.binding.basis.complexTypeDefinition):
    """Azimuth relative times for each internal calibration packet"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 987, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Element val uses Python identifier val
    __val = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'val'), 'val', '__aresysTypes_CTD_ANON_21_val', True, pyxb.utils.utility.Location('aresysTypes.xsd', 989, 16), )

    
    val = property(__val.value, __val.set, None, None)

    _ElementMap.update({
        __val.name() : __val
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.CTD_ANON_21 = CTD_ANON_21


# Complex type [anonymous] with content type ELEMENT_ONLY
class CTD_ANON_22 (pyxb.binding.basis.complexTypeDefinition):
    """Relative azimuth times for each SWL change"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1012, 14)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Element val uses Python identifier val
    __val = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'val'), 'val', '__aresysTypes_CTD_ANON_22_val', True, pyxb.utils.utility.Location('aresysTypes.xsd', 1014, 18), )

    
    val = property(__val.value, __val.set, None, None)

    _ElementMap.update({
        __val.name() : __val
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.CTD_ANON_22 = CTD_ANON_22


# Complex type [anonymous] with content type ELEMENT_ONLY
class CTD_ANON_23 (pyxb.binding.basis.complexTypeDefinition):
    """SWL changes values"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1031, 14)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Element val uses Python identifier val
    __val = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'val'), 'val', '__aresysTypes_CTD_ANON_23_val', True, pyxb.utils.utility.Location('aresysTypes.xsd', 1033, 18), )

    
    val = property(__val.value, __val.set, None, None)

    _ElementMap.update({
        __val.name() : __val
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.CTD_ANON_23 = CTD_ANON_23


# Complex type [anonymous] with content type ELEMENT_ONLY
class CTD_ANON_24 (pyxb.binding.basis.complexTypeDefinition):
    """Complex type [anonymous] with content type ELEMENT_ONLY"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1327, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Element DataBlockStatistic uses Python identifier DataBlockStatistic
    __DataBlockStatistic = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'DataBlockStatistic'), 'DataBlockStatistic', '__aresysTypes_CTD_ANON_24_DataBlockStatistic', True, pyxb.utils.utility.Location('aresysTypes.xsd', 1329, 16), )

    
    DataBlockStatistic = property(__DataBlockStatistic.value, __DataBlockStatistic.set, None, None)

    _ElementMap.update({
        __DataBlockStatistic.name() : __DataBlockStatistic
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.CTD_ANON_24 = CTD_ANON_24


# Complex type [anonymous] with content type ELEMENT_ONLY
class CTD_ANON_25 (pyxb.binding.basis.complexTypeDefinition):
    """Polynomial coefficients: const, rg, az, az*rg, rg^2, rg^3, rg^4 [Optional: rg^5 rg^6 .... rg^N]"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1469, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Element val uses Python identifier val
    __val = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'val'), 'val', '__aresysTypes_CTD_ANON_25_val', True, pyxb.utils.utility.Location('aresysTypes.xsd', 1471, 16), )

    
    val = property(__val.value, __val.set, None, None)

    _ElementMap.update({
        __val.name() : __val
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.CTD_ANON_25 = CTD_ANON_25


# Complex type [anonymous] with content type ELEMENT_ONLY
class CTD_ANON_26 (pyxb.binding.basis.complexTypeDefinition):
    """Polynomial coefficients of Range"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1523, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Element val uses Python identifier val
    __val = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'val'), 'val', '__aresysTypes_CTD_ANON_26_val', True, pyxb.utils.utility.Location('aresysTypes.xsd', 1525, 16), )

    
    val = property(__val.value, __val.set, None, None)

    _ElementMap.update({
        __val.name() : __val
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.CTD_ANON_26 = CTD_ANON_26


# Complex type [anonymous] with content type ELEMENT_ONLY
class CTD_ANON_27 (pyxb.binding.basis.complexTypeDefinition):
    """Polynomial coefficients of Azimuth"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1542, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Element val uses Python identifier val
    __val = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'val'), 'val', '__aresysTypes_CTD_ANON_27_val', True, pyxb.utils.utility.Location('aresysTypes.xsd', 1544, 16), )

    
    val = property(__val.value, __val.set, None, None)

    _ElementMap.update({
        __val.name() : __val
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.CTD_ANON_27 = CTD_ANON_27


# Complex type [anonymous] with content type ELEMENT_ONLY
class CTD_ANON_28 (pyxb.binding.basis.complexTypeDefinition):
    """Complex type [anonymous] with content type ELEMENT_ONLY"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1608, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Element Point uses Python identifier Point
    __Point = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'Point'), 'Point', '__aresysTypes_CTD_ANON_28_Point', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1610, 16), )

    
    Point = property(__Point.value, __Point.set, None, None)

    _ElementMap.update({
        __Point.name() : __Point
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.CTD_ANON_28 = CTD_ANON_28


# Complex type [anonymous] with content type ELEMENT_ONLY
class CTD_ANON_29 (pyxb.binding.basis.complexTypeDefinition):
    """Complex type [anonymous] with content type ELEMENT_ONLY"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1615, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Element Point uses Python identifier Point
    __Point = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'Point'), 'Point', '__aresysTypes_CTD_ANON_29_Point', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1617, 16), )

    
    Point = property(__Point.value, __Point.set, None, None)

    _ElementMap.update({
        __Point.name() : __Point
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.CTD_ANON_29 = CTD_ANON_29


# Complex type [anonymous] with content type ELEMENT_ONLY
class CTD_ANON_30 (pyxb.binding.basis.complexTypeDefinition):
    """Complex type [anonymous] with content type ELEMENT_ONLY"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1622, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Element Point uses Python identifier Point
    __Point = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'Point'), 'Point', '__aresysTypes_CTD_ANON_30_Point', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1624, 16), )

    
    Point = property(__Point.value, __Point.set, None, None)

    _ElementMap.update({
        __Point.name() : __Point
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.CTD_ANON_30 = CTD_ANON_30


# Complex type [anonymous] with content type ELEMENT_ONLY
class CTD_ANON_31 (pyxb.binding.basis.complexTypeDefinition):
    """Complex type [anonymous] with content type ELEMENT_ONLY"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1629, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Element Point uses Python identifier Point
    __Point = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'Point'), 'Point', '__aresysTypes_CTD_ANON_31_Point', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1631, 16), )

    
    Point = property(__Point.value, __Point.set, None, None)

    _ElementMap.update({
        __Point.name() : __Point
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.CTD_ANON_31 = CTD_ANON_31


# Complex type [anonymous] with content type ELEMENT_ONLY
class CTD_ANON_32 (pyxb.binding.basis.complexTypeDefinition):
    """Complex type [anonymous] with content type ELEMENT_ONLY"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1636, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Element Point uses Python identifier Point
    __Point = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'Point'), 'Point', '__aresysTypes_CTD_ANON_32_Point', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1638, 16), )

    
    Point = property(__Point.value, __Point.set, None, None)

    _ElementMap.update({
        __Point.name() : __Point
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.CTD_ANON_32 = CTD_ANON_32


# Complex type {aresysTypes}PointType with content type ELEMENT_ONLY
class PointType (pyxb.binding.basis.complexTypeDefinition):
    """Complex type {aresysTypes}PointType with content type ELEMENT_ONLY"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'PointType')
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1646, 2)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Element val uses Python identifier val
    __val = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'val'), 'val', '__aresysTypes_PointType_val', True, pyxb.utils.utility.Location('aresysTypes.xsd', 1648, 6), )

    
    val = property(__val.value, __val.set, None, None)

    _ElementMap.update({
        __val.name() : __val
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.PointType = PointType
Namespace.addCategoryObject('typeBinding', 'PointType', PointType)


# Complex type [anonymous] with content type ELEMENT_ONLY
class CTD_ANON_33 (pyxb.binding.basis.complexTypeDefinition):
    """Complex type [anonymous] with content type ELEMENT_ONLY"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1678, 16)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Element Lines uses Python identifier Lines
    __Lines = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'Lines'), 'Lines', '__aresysTypes_CTD_ANON_33_Lines', True, pyxb.utils.utility.Location('aresysTypes.xsd', 1680, 24), )

    
    Lines = property(__Lines.value, __Lines.set, None, None)

    _ElementMap.update({
        __Lines.name() : __Lines
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.CTD_ANON_33 = CTD_ANON_33


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_34 (pyxb.binding.basis.complexTypeDefinition):
    """Complex type [anonymous] with content type SIMPLE"""
    _TypeDefinition = pyxb.binding.datatypes.unsignedInt
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1681, 28)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.unsignedInt
    
    # Attribute FromBurst uses Python identifier FromBurst
    __FromBurst = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'FromBurst'), 'FromBurst', '__aresysTypes_CTD_ANON_34_FromBurst', pyxb.binding.datatypes.unsignedInt, required=True)
    __FromBurst._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1684, 40)
    __FromBurst._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1684, 40)
    
    FromBurst = property(__FromBurst.value, __FromBurst.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __FromBurst.name() : __FromBurst
    })
_module_typeBindings.CTD_ANON_34 = CTD_ANON_34


# Complex type {aresysTypes}BurstType with content type ELEMENT_ONLY
class BurstType (pyxb.binding.basis.complexTypeDefinition):
    """Complex type {aresysTypes}BurstType with content type ELEMENT_ONLY"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'BurstType')
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1743, 2)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Element RangeStartTime uses Python identifier RangeStartTime
    __RangeStartTime = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'RangeStartTime'), 'RangeStartTime', '__aresysTypes_BurstType_RangeStartTime', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1745, 6), )

    
    RangeStartTime = property(__RangeStartTime.value, __RangeStartTime.set, None, 'Range absolute start time [s]')

    
    # Element AzimuthStartTime uses Python identifier AzimuthStartTime
    __AzimuthStartTime = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'AzimuthStartTime'), 'AzimuthStartTime', '__aresysTypes_BurstType_AzimuthStartTime', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1750, 6), )

    
    AzimuthStartTime = property(__AzimuthStartTime.value, __AzimuthStartTime.set, None, 'Azimuth start time absolute value [Utc]')

    
    # Element BurstCenterAzimuthShift uses Python identifier BurstCenterAzimuthShift
    __BurstCenterAzimuthShift = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'BurstCenterAzimuthShift'), 'BurstCenterAzimuthShift', '__aresysTypes_BurstType_BurstCenterAzimuthShift', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1755, 6), )

    
    BurstCenterAzimuthShift = property(__BurstCenterAzimuthShift.value, __BurstCenterAzimuthShift.set, None, None)

    
    # Attribute N uses Python identifier N
    __N = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'N'), 'N', '__aresysTypes_BurstType_N', pyxb.binding.datatypes.unsignedInt, required=True)
    __N._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1757, 4)
    __N._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1757, 4)
    
    N = property(__N.value, __N.set, None, None)

    _ElementMap.update({
        __RangeStartTime.name() : __RangeStartTime,
        __AzimuthStartTime.name() : __AzimuthStartTime,
        __BurstCenterAzimuthShift.name() : __BurstCenterAzimuthShift
    })
    _AttributeMap.update({
        __N.name() : __N
    })
_module_typeBindings.BurstType = BurstType
Namespace.addCategoryObject('typeBinding', 'BurstType', BurstType)


# Complex type {aresysTypes}DCOMPLEX with content type ELEMENT_ONLY
class DCOMPLEX (pyxb.binding.basis.complexTypeDefinition):
    """Complex type {aresysTypes}DCOMPLEX with content type ELEMENT_ONLY"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'DCOMPLEX')
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1759, 2)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Element RealValue uses Python identifier RealValue
    __RealValue = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'RealValue'), 'RealValue', '__aresysTypes_DCOMPLEX_RealValue', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1761, 6), )

    
    RealValue = property(__RealValue.value, __RealValue.set, None, None)

    
    # Element ImaginaryValue uses Python identifier ImaginaryValue
    __ImaginaryValue = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'ImaginaryValue'), 'ImaginaryValue', '__aresysTypes_DCOMPLEX_ImaginaryValue', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1762, 6), )

    
    ImaginaryValue = property(__ImaginaryValue.value, __ImaginaryValue.set, None, None)

    _ElementMap.update({
        __RealValue.name() : __RealValue,
        __ImaginaryValue.name() : __ImaginaryValue
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.DCOMPLEX = DCOMPLEX
Namespace.addCategoryObject('typeBinding', 'DCOMPLEX', DCOMPLEX)


# Complex type {aresysTypes}doubleWithUnit with content type SIMPLE
class doubleWithUnit (pyxb.binding.basis.complexTypeDefinition):
    """Complex type {aresysTypes}doubleWithUnit with content type SIMPLE"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'doubleWithUnit')
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 160, 2)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_doubleWithUnit_unit', _module_typeBindings.units, required=True)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 163, 8)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 163, 8)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.doubleWithUnit = doubleWithUnit
Namespace.addCategoryObject('typeBinding', 'doubleWithUnit', doubleWithUnit)


# Complex type {aresysTypes}stringWithUnit with content type SIMPLE
class stringWithUnit (pyxb.binding.basis.complexTypeDefinition):
    """Complex type {aresysTypes}stringWithUnit with content type SIMPLE"""
    _TypeDefinition = pyxb.binding.datatypes.string
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'stringWithUnit')
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 167, 2)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.string
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_stringWithUnit_unit', _module_typeBindings.units, required=True)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 170, 8)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 170, 8)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.stringWithUnit = stringWithUnit
Namespace.addCategoryObject('typeBinding', 'stringWithUnit', stringWithUnit)


# Complex type {aresysTypes}ROIType with content type ELEMENT_ONLY
class ROIType (TreeElementBaseType):
    """Information regarding the region of interest related to the raster data"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'ROIType')
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 174, 2)
    _ElementMap = TreeElementBaseType._ElementMap.copy()
    _AttributeMap = TreeElementBaseType._AttributeMap.copy()
    # Base type is TreeElementBaseType
    
    # Element SlowStartTime uses Python identifier SlowStartTime
    __SlowStartTime = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'SlowStartTime'), 'SlowStartTime', '__aresysTypes_ROIType_SlowStartTime', False, pyxb.utils.utility.Location('aresysTypes.xsd', 181, 10), )

    
    SlowStartTime = property(__SlowStartTime.value, __SlowStartTime.set, None, None)

    
    # Element SlowTimeLength uses Python identifier SlowTimeLength
    __SlowTimeLength = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'SlowTimeLength'), 'SlowTimeLength', '__aresysTypes_ROIType_SlowTimeLength', False, pyxb.utils.utility.Location('aresysTypes.xsd', 182, 10), )

    
    SlowTimeLength = property(__SlowTimeLength.value, __SlowTimeLength.set, None, None)

    
    # Element SlowTimeDelta uses Python identifier SlowTimeDelta
    __SlowTimeDelta = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'SlowTimeDelta'), 'SlowTimeDelta', '__aresysTypes_ROIType_SlowTimeDelta', False, pyxb.utils.utility.Location('aresysTypes.xsd', 183, 10), )

    
    SlowTimeDelta = property(__SlowTimeDelta.value, __SlowTimeDelta.set, None, None)

    
    # Element FastStartTime uses Python identifier FastStartTime
    __FastStartTime = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'FastStartTime'), 'FastStartTime', '__aresysTypes_ROIType_FastStartTime', False, pyxb.utils.utility.Location('aresysTypes.xsd', 184, 10), )

    
    FastStartTime = property(__FastStartTime.value, __FastStartTime.set, None, None)

    
    # Element FastTimeLength uses Python identifier FastTimeLength
    __FastTimeLength = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'FastTimeLength'), 'FastTimeLength', '__aresysTypes_ROIType_FastTimeLength', False, pyxb.utils.utility.Location('aresysTypes.xsd', 191, 10), )

    
    FastTimeLength = property(__FastTimeLength.value, __FastTimeLength.set, None, None)

    
    # Element FastTimeDelta uses Python identifier FastTimeDelta
    __FastTimeDelta = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'FastTimeDelta'), 'FastTimeDelta', '__aresysTypes_ROIType_FastTimeDelta', False, pyxb.utils.utility.Location('aresysTypes.xsd', 198, 10), )

    
    FastTimeDelta = property(__FastTimeDelta.value, __FastTimeDelta.set, None, None)

    
    # Element Lines uses Python identifier Lines
    __Lines = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'Lines'), 'Lines', '__aresysTypes_ROIType_Lines', False, pyxb.utils.utility.Location('aresysTypes.xsd', 205, 10), )

    
    Lines = property(__Lines.value, __Lines.set, None, None)

    
    # Element Samples uses Python identifier Samples
    __Samples = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'Samples'), 'Samples', '__aresysTypes_ROIType_Samples', False, pyxb.utils.utility.Location('aresysTypes.xsd', 212, 10), )

    
    Samples = property(__Samples.value, __Samples.set, None, None)

    
    # Attribute Number inherited from {aresysTypes}TreeElementBaseType
    
    # Attribute Total inherited from {aresysTypes}TreeElementBaseType
    _ElementMap.update({
        __SlowStartTime.name() : __SlowStartTime,
        __SlowTimeLength.name() : __SlowTimeLength,
        __SlowTimeDelta.name() : __SlowTimeDelta,
        __FastStartTime.name() : __FastStartTime,
        __FastTimeLength.name() : __FastTimeLength,
        __FastTimeDelta.name() : __FastTimeDelta,
        __Lines.name() : __Lines,
        __Samples.name() : __Samples
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.ROIType = ROIType
Namespace.addCategoryObject('typeBinding', 'ROIType', ROIType)


# Complex type {aresysTypes}PulseType with content type ELEMENT_ONLY
class PulseType (TreeElementBaseType):
    """Transmitted pulse parameters"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'PulseType')
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 223, 2)
    _ElementMap = TreeElementBaseType._ElementMap.copy()
    _AttributeMap = TreeElementBaseType._AttributeMap.copy()
    # Base type is TreeElementBaseType
    
    # Element Direction uses Python identifier Direction
    __Direction = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'Direction'), 'Direction', '__aresysTypes_PulseType_Direction', False, pyxb.utils.utility.Location('aresysTypes.xsd', 230, 10), )

    
    Direction = property(__Direction.value, __Direction.set, None, 'Pulse direction (UP, DOWN)')

    
    # Element PulseLength uses Python identifier PulseLength
    __PulseLength = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'PulseLength'), 'PulseLength', '__aresysTypes_PulseType_PulseLength', False, pyxb.utils.utility.Location('aresysTypes.xsd', 241, 10), )

    
    PulseLength = property(__PulseLength.value, __PulseLength.set, None, 'Pulse length [s]')

    
    # Element Bandwidth uses Python identifier Bandwidth
    __Bandwidth = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'Bandwidth'), 'Bandwidth', '__aresysTypes_PulseType_Bandwidth', False, pyxb.utils.utility.Location('aresysTypes.xsd', 246, 10), )

    
    Bandwidth = property(__Bandwidth.value, __Bandwidth.set, None, 'Pulse bandwidth [Hz]')

    
    # Element PulseEnergy uses Python identifier PulseEnergy
    __PulseEnergy = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'PulseEnergy'), 'PulseEnergy', '__aresysTypes_PulseType_PulseEnergy', False, pyxb.utils.utility.Location('aresysTypes.xsd', 251, 10), )

    
    PulseEnergy = property(__PulseEnergy.value, __PulseEnergy.set, None, 'Pulse energy [J]')

    
    # Element PulseSamplingRate uses Python identifier PulseSamplingRate
    __PulseSamplingRate = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'PulseSamplingRate'), 'PulseSamplingRate', '__aresysTypes_PulseType_PulseSamplingRate', False, pyxb.utils.utility.Location('aresysTypes.xsd', 256, 10), )

    
    PulseSamplingRate = property(__PulseSamplingRate.value, __PulseSamplingRate.set, None, 'Pulse sampling rate [Hz]')

    
    # Element PulseStartFrequency uses Python identifier PulseStartFrequency
    __PulseStartFrequency = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'PulseStartFrequency'), 'PulseStartFrequency', '__aresysTypes_PulseType_PulseStartFrequency', False, pyxb.utils.utility.Location('aresysTypes.xsd', 261, 10), )

    
    PulseStartFrequency = property(__PulseStartFrequency.value, __PulseStartFrequency.set, None, 'Pulse start frequency [Hz]')

    
    # Element PulseStartPhase uses Python identifier PulseStartPhase
    __PulseStartPhase = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'PulseStartPhase'), 'PulseStartPhase', '__aresysTypes_PulseType_PulseStartPhase', False, pyxb.utils.utility.Location('aresysTypes.xsd', 266, 10), )

    
    PulseStartPhase = property(__PulseStartPhase.value, __PulseStartPhase.set, None, 'Pulse start phase [rad]')

    
    # Attribute Number inherited from {aresysTypes}TreeElementBaseType
    
    # Attribute Total inherited from {aresysTypes}TreeElementBaseType
    _ElementMap.update({
        __Direction.name() : __Direction,
        __PulseLength.name() : __PulseLength,
        __Bandwidth.name() : __Bandwidth,
        __PulseEnergy.name() : __PulseEnergy,
        __PulseSamplingRate.name() : __PulseSamplingRate,
        __PulseStartFrequency.name() : __PulseStartFrequency,
        __PulseStartPhase.name() : __PulseStartPhase
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.PulseType = PulseType
Namespace.addCategoryObject('typeBinding', 'PulseType', PulseType)


# Complex type {aresysTypes}DataSetInfoType with content type ELEMENT_ONLY
class DataSetInfoType (TreeElementBaseType):
    """Information regarding the dataset"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'DataSetInfoType')
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 275, 2)
    _ElementMap = TreeElementBaseType._ElementMap.copy()
    _AttributeMap = TreeElementBaseType._AttributeMap.copy()
    # Base type is TreeElementBaseType
    
    # Element SensorName uses Python identifier SensorName
    __SensorName = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'SensorName'), 'SensorName', '__aresysTypes_DataSetInfoType_SensorName', False, pyxb.utils.utility.Location('aresysTypes.xsd', 282, 10), )

    
    SensorName = property(__SensorName.value, __SensorName.set, None, 'Name of the sensor used to acquire the image: ASAR, PALSAR, ...')

    
    # Element Description uses Python identifier Description
    __Description = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'Description'), 'Description', '__aresysTypes_DataSetInfoType_Description', False, pyxb.utils.utility.Location('aresysTypes.xsd', 287, 10), )

    
    Description = property(__Description.value, __Description.set, None, 'Description of the image')

    
    # Element SenseDate uses Python identifier SenseDate
    __SenseDate = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'SenseDate'), 'SenseDate', '__aresysTypes_DataSetInfoType_SenseDate', False, pyxb.utils.utility.Location('aresysTypes.xsd', 299, 10), )

    
    SenseDate = property(__SenseDate.value, __SenseDate.set, None, 'Image acquisition date')

    
    # Element AcquisitionMode uses Python identifier AcquisitionMode
    __AcquisitionMode = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'AcquisitionMode'), 'AcquisitionMode', '__aresysTypes_DataSetInfoType_AcquisitionMode', False, pyxb.utils.utility.Location('aresysTypes.xsd', 311, 10), )

    
    AcquisitionMode = property(__AcquisitionMode.value, __AcquisitionMode.set, None, 'Image acquisition mode: STRIPMAP, TOPSAR, ...')

    
    # Element ImageType uses Python identifier ImageType
    __ImageType = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'ImageType'), 'ImageType', '__aresysTypes_DataSetInfoType_ImageType', False, pyxb.utils.utility.Location('aresysTypes.xsd', 323, 10), )

    
    ImageType = property(__ImageType.value, __ImageType.set, None, 'Image type: RAW DATA, RANGE FOCUSED, AZIMUTH FOCUSED ')

    
    # Element Projection uses Python identifier Projection
    __Projection = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'Projection'), 'Projection', '__aresysTypes_DataSetInfoType_Projection', False, pyxb.utils.utility.Location('aresysTypes.xsd', 335, 10), )

    
    Projection = property(__Projection.value, __Projection.set, None, 'Image projection: SLANT RANGE, GROUND RANGE')

    
    # Element ProjectionParameters uses Python identifier ProjectionParameters
    __ProjectionParameters = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'ProjectionParameters'), 'ProjectionParameters', '__aresysTypes_DataSetInfoType_ProjectionParameters', False, pyxb.utils.utility.Location('aresysTypes.xsd', 347, 10), )

    
    ProjectionParameters = property(__ProjectionParameters.value, __ProjectionParameters.set, None, None)

    
    # Element AcquisitionStation uses Python identifier AcquisitionStation
    __AcquisitionStation = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'AcquisitionStation'), 'AcquisitionStation', '__aresysTypes_DataSetInfoType_AcquisitionStation', False, pyxb.utils.utility.Location('aresysTypes.xsd', 356, 10), )

    
    AcquisitionStation = property(__AcquisitionStation.value, __AcquisitionStation.set, None, 'Image acquisition station')

    
    # Element ProcessingCenter uses Python identifier ProcessingCenter
    __ProcessingCenter = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'ProcessingCenter'), 'ProcessingCenter', '__aresysTypes_DataSetInfoType_ProcessingCenter', False, pyxb.utils.utility.Location('aresysTypes.xsd', 368, 10), )

    
    ProcessingCenter = property(__ProcessingCenter.value, __ProcessingCenter.set, None, 'Image processing center')

    
    # Element ProcessingDate uses Python identifier ProcessingDate
    __ProcessingDate = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'ProcessingDate'), 'ProcessingDate', '__aresysTypes_DataSetInfoType_ProcessingDate', False, pyxb.utils.utility.Location('aresysTypes.xsd', 380, 10), )

    
    ProcessingDate = property(__ProcessingDate.value, __ProcessingDate.set, None, 'Image processing date')

    
    # Element ProcessingSoftware uses Python identifier ProcessingSoftware
    __ProcessingSoftware = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'ProcessingSoftware'), 'ProcessingSoftware', '__aresysTypes_DataSetInfoType_ProcessingSoftware', False, pyxb.utils.utility.Location('aresysTypes.xsd', 392, 10), )

    
    ProcessingSoftware = property(__ProcessingSoftware.value, __ProcessingSoftware.set, None, 'Image processing software')

    
    # Element fc_hz uses Python identifier fc_hz
    __fc_hz = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'fc_hz'), 'fc_hz', '__aresysTypes_DataSetInfoType_fc_hz', False, pyxb.utils.utility.Location('aresysTypes.xsd', 404, 10), )

    
    fc_hz = property(__fc_hz.value, __fc_hz.set, None, 'Radar carrier frequency [Hz]')

    
    # Element SideLooking uses Python identifier SideLooking
    __SideLooking = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'SideLooking'), 'SideLooking', '__aresysTypes_DataSetInfoType_SideLooking', False, pyxb.utils.utility.Location('aresysTypes.xsd', 416, 10), )

    
    SideLooking = property(__SideLooking.value, __SideLooking.set, None, 'Radar side looking: LEFT, RIGHT')

    
    # Element ExternalCalibrationFactor uses Python identifier ExternalCalibrationFactor
    __ExternalCalibrationFactor = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'ExternalCalibrationFactor'), 'ExternalCalibrationFactor', '__aresysTypes_DataSetInfoType_ExternalCalibrationFactor', False, pyxb.utils.utility.Location('aresysTypes.xsd', 421, 10), )

    
    ExternalCalibrationFactor = property(__ExternalCalibrationFactor.value, __ExternalCalibrationFactor.set, None, 'External calibration factor')

    
    # Element DataTakeID uses Python identifier DataTakeID
    __DataTakeID = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'DataTakeID'), 'DataTakeID', '__aresysTypes_DataSetInfoType_DataTakeID', False, pyxb.utils.utility.Location('aresysTypes.xsd', 426, 10), )

    
    DataTakeID = property(__DataTakeID.value, __DataTakeID.set, None, None)

    
    # Element InstrumentConfID uses Python identifier InstrumentConfID
    __InstrumentConfID = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'InstrumentConfID'), 'InstrumentConfID', '__aresysTypes_DataSetInfoType_InstrumentConfID', False, pyxb.utils.utility.Location('aresysTypes.xsd', 427, 10), )

    
    InstrumentConfID = property(__InstrumentConfID.value, __InstrumentConfID.set, None, None)

    
    # Attribute Number inherited from {aresysTypes}TreeElementBaseType
    
    # Attribute Total inherited from {aresysTypes}TreeElementBaseType
    _ElementMap.update({
        __SensorName.name() : __SensorName,
        __Description.name() : __Description,
        __SenseDate.name() : __SenseDate,
        __AcquisitionMode.name() : __AcquisitionMode,
        __ImageType.name() : __ImageType,
        __Projection.name() : __Projection,
        __ProjectionParameters.name() : __ProjectionParameters,
        __AcquisitionStation.name() : __AcquisitionStation,
        __ProcessingCenter.name() : __ProcessingCenter,
        __ProcessingDate.name() : __ProcessingDate,
        __ProcessingSoftware.name() : __ProcessingSoftware,
        __fc_hz.name() : __fc_hz,
        __SideLooking.name() : __SideLooking,
        __ExternalCalibrationFactor.name() : __ExternalCalibrationFactor,
        __DataTakeID.name() : __DataTakeID,
        __InstrumentConfID.name() : __InstrumentConfID
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.DataSetInfoType = DataSetInfoType
Namespace.addCategoryObject('typeBinding', 'DataSetInfoType', DataSetInfoType)


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_35 (pyxb.binding.basis.complexTypeDefinition):
    """Description of the image"""
    _TypeDefinition = pyxb.binding.datatypes.string
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 291, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.string
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_35_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 294, 18)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 294, 18)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_35 = CTD_ANON_35


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_36 (pyxb.binding.basis.complexTypeDefinition):
    """Image acquisition date"""
    _TypeDefinition = pyxb.binding.datatypes.string
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 303, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.string
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_36_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 306, 18)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 306, 18)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_36 = CTD_ANON_36


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_37 (pyxb.binding.basis.complexTypeDefinition):
    """Image acquisition mode: STRIPMAP, TOPSAR, ..."""
    _TypeDefinition = pyxb.binding.datatypes.string
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 315, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.string
    
    # Attribute Polarization uses Python identifier Polarization
    __Polarization = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'Polarization'), 'Polarization', '__aresysTypes_CTD_ANON_37_Polarization', _module_typeBindings.GlobalPolarizationType)
    __Polarization._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 318, 18)
    __Polarization._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 318, 18)
    
    Polarization = property(__Polarization.value, __Polarization.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __Polarization.name() : __Polarization
    })
_module_typeBindings.CTD_ANON_37 = CTD_ANON_37


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_38 (pyxb.binding.basis.complexTypeDefinition):
    """Image type: RAW DATA, RANGE FOCUSED, AZIMUTH FOCUSED """
    _TypeDefinition = pyxb.binding.datatypes.string
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 327, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.string
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_38_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 330, 18)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 330, 18)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_38 = CTD_ANON_38


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_39 (pyxb.binding.basis.complexTypeDefinition):
    """Image projection: SLANT RANGE, GROUND RANGE"""
    _TypeDefinition = pyxb.binding.datatypes.string
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 339, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.string
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_39_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 342, 18)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 342, 18)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_39 = CTD_ANON_39


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_40 (pyxb.binding.basis.complexTypeDefinition):
    """Image acquisition station"""
    _TypeDefinition = pyxb.binding.datatypes.string
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 360, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.string
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_40_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 363, 18)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 363, 18)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_40 = CTD_ANON_40


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_41 (pyxb.binding.basis.complexTypeDefinition):
    """Image processing center"""
    _TypeDefinition = pyxb.binding.datatypes.string
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 372, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.string
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_41_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 375, 18)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 375, 18)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_41 = CTD_ANON_41


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_42 (pyxb.binding.basis.complexTypeDefinition):
    """Image processing date"""
    _TypeDefinition = pyxb.binding.datatypes.string
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 384, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.string
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_42_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 387, 18)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 387, 18)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_42 = CTD_ANON_42


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_43 (pyxb.binding.basis.complexTypeDefinition):
    """Image processing software"""
    _TypeDefinition = pyxb.binding.datatypes.string
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 396, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.string
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_43_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 399, 18)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 399, 18)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_43 = CTD_ANON_43


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_44 (pyxb.binding.basis.complexTypeDefinition):
    """Radar carrier frequency [Hz]"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 408, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_44_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 411, 18)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 411, 18)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_44 = CTD_ANON_44


# Complex type {aresysTypes}StateVectorDataType with content type ELEMENT_ONLY
class StateVectorDataType (TreeElementBaseType):
    """Information regarding position and velocity of the sensor along the orbit (State Vector Data)"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'StateVectorDataType')
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 439, 2)
    _ElementMap = TreeElementBaseType._ElementMap.copy()
    _AttributeMap = TreeElementBaseType._AttributeMap.copy()
    # Base type is TreeElementBaseType
    
    # Element OrbitNumber uses Python identifier OrbitNumber
    __OrbitNumber = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'OrbitNumber'), 'OrbitNumber', '__aresysTypes_StateVectorDataType_OrbitNumber', False, pyxb.utils.utility.Location('aresysTypes.xsd', 446, 10), )

    
    OrbitNumber = property(__OrbitNumber.value, __OrbitNumber.set, None, 'Number of the orbit')

    
    # Element Track uses Python identifier Track
    __Track = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'Track'), 'Track', '__aresysTypes_StateVectorDataType_Track', False, pyxb.utils.utility.Location('aresysTypes.xsd', 451, 10), )

    
    Track = property(__Track.value, __Track.set, None, 'Number of the track')

    
    # Element OrbitDirection uses Python identifier OrbitDirection
    __OrbitDirection = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'OrbitDirection'), 'OrbitDirection', '__aresysTypes_StateVectorDataType_OrbitDirection', False, pyxb.utils.utility.Location('aresysTypes.xsd', 456, 10), )

    
    OrbitDirection = property(__OrbitDirection.value, __OrbitDirection.set, None, 'Direction of the orbit: ASCENDING, DESCENDING')

    
    # Element pSV_m uses Python identifier pSV_m
    __pSV_m = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'pSV_m'), 'pSV_m', '__aresysTypes_StateVectorDataType_pSV_m', False, pyxb.utils.utility.Location('aresysTypes.xsd', 466, 10), )

    
    pSV_m = property(__pSV_m.value, __pSV_m.set, None, 'Orbit state vectors position coordinates (xyz) [m]')

    
    # Element vSV_mOs uses Python identifier vSV_mOs
    __vSV_mOs = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'vSV_mOs'), 'vSV_mOs', '__aresysTypes_StateVectorDataType_vSV_mOs', False, pyxb.utils.utility.Location('aresysTypes.xsd', 485, 10), )

    
    vSV_mOs = property(__vSV_mOs.value, __vSV_mOs.set, None, 'Orbit state vectors velocity coordinates [m/s]')

    
    # Element t_ref_Utc uses Python identifier t_ref_Utc
    __t_ref_Utc = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 't_ref_Utc'), 't_ref_Utc', '__aresysTypes_StateVectorDataType_t_ref_Utc', False, pyxb.utils.utility.Location('aresysTypes.xsd', 504, 10), )

    
    t_ref_Utc = property(__t_ref_Utc.value, __t_ref_Utc.set, None, 'Azimuth absolute start time for the first state vector [Utc]')

    
    # Element dtSV_s uses Python identifier dtSV_s
    __dtSV_s = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'dtSV_s'), 'dtSV_s', '__aresysTypes_StateVectorDataType_dtSV_s', False, pyxb.utils.utility.Location('aresysTypes.xsd', 509, 10), )

    
    dtSV_s = property(__dtSV_s.value, __dtSV_s.set, None, 'Azimuth time interval between two consecutive state vectors [s]')

    
    # Element nSV_n uses Python identifier nSV_n
    __nSV_n = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'nSV_n'), 'nSV_n', '__aresysTypes_StateVectorDataType_nSV_n', False, pyxb.utils.utility.Location('aresysTypes.xsd', 521, 10), )

    
    nSV_n = property(__nSV_n.value, __nSV_n.set, None, 'Number of state vectors')

    
    # Element AscendingNodeTime uses Python identifier AscendingNodeTime
    __AscendingNodeTime = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'AscendingNodeTime'), 'AscendingNodeTime', '__aresysTypes_StateVectorDataType_AscendingNodeTime', False, pyxb.utils.utility.Location('aresysTypes.xsd', 533, 10), )

    
    AscendingNodeTime = property(__AscendingNodeTime.value, __AscendingNodeTime.set, None, 'Azimuth absolute time of the ascending node')

    
    # Element AscendingNodeCoords uses Python identifier AscendingNodeCoords
    __AscendingNodeCoords = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'AscendingNodeCoords'), 'AscendingNodeCoords', '__aresysTypes_StateVectorDataType_AscendingNodeCoords', False, pyxb.utils.utility.Location('aresysTypes.xsd', 538, 10), )

    
    AscendingNodeCoords = property(__AscendingNodeCoords.value, __AscendingNodeCoords.set, None, 'Coordinates of the ascending node')

    
    # Attribute Number inherited from {aresysTypes}TreeElementBaseType
    
    # Attribute Total inherited from {aresysTypes}TreeElementBaseType
    _ElementMap.update({
        __OrbitNumber.name() : __OrbitNumber,
        __Track.name() : __Track,
        __OrbitDirection.name() : __OrbitDirection,
        __pSV_m.name() : __pSV_m,
        __vSV_mOs.name() : __vSV_mOs,
        __t_ref_Utc.name() : __t_ref_Utc,
        __dtSV_s.name() : __dtSV_s,
        __nSV_n.name() : __nSV_n,
        __AscendingNodeTime.name() : __AscendingNodeTime,
        __AscendingNodeCoords.name() : __AscendingNodeCoords
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.StateVectorDataType = StateVectorDataType
Namespace.addCategoryObject('typeBinding', 'StateVectorDataType', StateVectorDataType)


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_45 (pyxb.binding.basis.complexTypeDefinition):
    """Direction of the orbit: ASCENDING, DESCENDING"""
    _TypeDefinition = AscendingDescendingType
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 460, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is AscendingDescendingType
    _ElementMap.update({
        
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.CTD_ANON_45 = CTD_ANON_45


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_46 (pyxb.binding.basis.complexTypeDefinition):
    """Complex type [anonymous] with content type SIMPLE"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 473, 18)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute N uses Python identifier N
    __N = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'N'), 'N', '__aresysTypes_CTD_ANON_46_N', pyxb.binding.datatypes.int)
    __N._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 476, 24)
    __N._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 476, 24)
    
    N = property(__N.value, __N.set, None, None)

    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_46_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 477, 24)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 477, 24)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __N.name() : __N,
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_46 = CTD_ANON_46


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_47 (pyxb.binding.basis.complexTypeDefinition):
    """Complex type [anonymous] with content type SIMPLE"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 492, 18)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute N uses Python identifier N
    __N = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'N'), 'N', '__aresysTypes_CTD_ANON_47_N', pyxb.binding.datatypes.int)
    __N._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 495, 24)
    __N._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 495, 24)
    
    N = property(__N.value, __N.set, None, None)

    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_47_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 496, 24)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 496, 24)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __N.name() : __N,
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_47 = CTD_ANON_47


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_48 (pyxb.binding.basis.complexTypeDefinition):
    """Azimuth time interval between two consecutive state vectors [s]"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 513, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_48_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 516, 18)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 516, 18)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_48 = CTD_ANON_48


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_49 (pyxb.binding.basis.complexTypeDefinition):
    """Number of state vectors"""
    _TypeDefinition = pyxb.binding.datatypes.int
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 525, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.int
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_49_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 528, 18)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 528, 18)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_49 = CTD_ANON_49


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_50 (pyxb.binding.basis.complexTypeDefinition):
    """Complex type [anonymous] with content type SIMPLE"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 545, 18)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute N uses Python identifier N
    __N = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'N'), 'N', '__aresysTypes_CTD_ANON_50_N', pyxb.binding.datatypes.int)
    __N._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 548, 24)
    __N._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 548, 24)
    
    N = property(__N.value, __N.set, None, None)

    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_50_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 549, 24)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 549, 24)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __N.name() : __N,
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_50 = CTD_ANON_50


# Complex type {aresysTypes}AttitudeInfoType with content type ELEMENT_ONLY
class AttitudeInfoType (TreeElementBaseType):
    """Sensor attitude information"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'AttitudeInfoType')
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 561, 2)
    _ElementMap = TreeElementBaseType._ElementMap.copy()
    _AttributeMap = TreeElementBaseType._AttributeMap.copy()
    # Base type is TreeElementBaseType
    
    # Element t_ref_Utc uses Python identifier t_ref_Utc
    __t_ref_Utc = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 't_ref_Utc'), 't_ref_Utc', '__aresysTypes_AttitudeInfoType_t_ref_Utc', False, pyxb.utils.utility.Location('aresysTypes.xsd', 568, 10), )

    
    t_ref_Utc = property(__t_ref_Utc.value, __t_ref_Utc.set, None, 'Azimuth absolute start time for the first attitude value [Utc]')

    
    # Element dtYPR_s uses Python identifier dtYPR_s
    __dtYPR_s = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'dtYPR_s'), 'dtYPR_s', '__aresysTypes_AttitudeInfoType_dtYPR_s', False, pyxb.utils.utility.Location('aresysTypes.xsd', 573, 10), )

    
    dtYPR_s = property(__dtYPR_s.value, __dtYPR_s.set, None, 'Azimuth time interval between two consecutive attitude values [s]')

    
    # Element nYPR_n uses Python identifier nYPR_n
    __nYPR_n = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'nYPR_n'), 'nYPR_n', '__aresysTypes_AttitudeInfoType_nYPR_n', False, pyxb.utils.utility.Location('aresysTypes.xsd', 585, 10), )

    
    nYPR_n = property(__nYPR_n.value, __nYPR_n.set, None, 'Number of attitude values')

    
    # Element yaw_deg uses Python identifier yaw_deg
    __yaw_deg = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'yaw_deg'), 'yaw_deg', '__aresysTypes_AttitudeInfoType_yaw_deg', False, pyxb.utils.utility.Location('aresysTypes.xsd', 597, 10), )

    
    yaw_deg = property(__yaw_deg.value, __yaw_deg.set, None, 'Yaw angle values [deg]')

    
    # Element pitch_deg uses Python identifier pitch_deg
    __pitch_deg = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'pitch_deg'), 'pitch_deg', '__aresysTypes_AttitudeInfoType_pitch_deg', False, pyxb.utils.utility.Location('aresysTypes.xsd', 616, 10), )

    
    pitch_deg = property(__pitch_deg.value, __pitch_deg.set, None, 'Pitch angle values [deg]')

    
    # Element roll_deg uses Python identifier roll_deg
    __roll_deg = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'roll_deg'), 'roll_deg', '__aresysTypes_AttitudeInfoType_roll_deg', False, pyxb.utils.utility.Location('aresysTypes.xsd', 635, 10), )

    
    roll_deg = property(__roll_deg.value, __roll_deg.set, None, 'Roll angle values [deg]')

    
    # Element referenceFrame uses Python identifier referenceFrame
    __referenceFrame = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'referenceFrame'), 'referenceFrame', '__aresysTypes_AttitudeInfoType_referenceFrame', False, pyxb.utils.utility.Location('aresysTypes.xsd', 654, 10), )

    
    referenceFrame = property(__referenceFrame.value, __referenceFrame.set, None, 'Reference frame: GEOCENTRIC, GEODETIC, ZERODOPPLER')

    
    # Element rotationOrder uses Python identifier rotationOrder
    __rotationOrder = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'rotationOrder'), 'rotationOrder', '__aresysTypes_AttitudeInfoType_rotationOrder', False, pyxb.utils.utility.Location('aresysTypes.xsd', 659, 10), )

    
    rotationOrder = property(__rotationOrder.value, __rotationOrder.set, None, 'Rotation order: YPR, YRP, PRY, PYR, RPY, RYP')

    
    # Element AttitudeType uses Python identifier AttitudeType
    __AttitudeType = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'AttitudeType'), 'AttitudeType', '__aresysTypes_AttitudeInfoType_AttitudeType', False, pyxb.utils.utility.Location('aresysTypes.xsd', 664, 10), )

    
    AttitudeType = property(__AttitudeType.value, __AttitudeType.set, None, 'Attitude type: NOMINAL, REFINED')

    
    # Attribute Number inherited from {aresysTypes}TreeElementBaseType
    
    # Attribute Total inherited from {aresysTypes}TreeElementBaseType
    _ElementMap.update({
        __t_ref_Utc.name() : __t_ref_Utc,
        __dtYPR_s.name() : __dtYPR_s,
        __nYPR_n.name() : __nYPR_n,
        __yaw_deg.name() : __yaw_deg,
        __pitch_deg.name() : __pitch_deg,
        __roll_deg.name() : __roll_deg,
        __referenceFrame.name() : __referenceFrame,
        __rotationOrder.name() : __rotationOrder,
        __AttitudeType.name() : __AttitudeType
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.AttitudeInfoType = AttitudeInfoType
Namespace.addCategoryObject('typeBinding', 'AttitudeInfoType', AttitudeInfoType)


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_51 (pyxb.binding.basis.complexTypeDefinition):
    """Azimuth time interval between two consecutive attitude values [s]"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 577, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_51_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 580, 18)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 580, 18)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_51 = CTD_ANON_51


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_52 (pyxb.binding.basis.complexTypeDefinition):
    """Number of attitude values"""
    _TypeDefinition = pyxb.binding.datatypes.int
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 589, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.int
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_52_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 592, 18)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 592, 18)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_52 = CTD_ANON_52


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_53 (pyxb.binding.basis.complexTypeDefinition):
    """Complex type [anonymous] with content type SIMPLE"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 604, 18)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute N uses Python identifier N
    __N = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'N'), 'N', '__aresysTypes_CTD_ANON_53_N', pyxb.binding.datatypes.int)
    __N._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 607, 24)
    __N._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 607, 24)
    
    N = property(__N.value, __N.set, None, None)

    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_53_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 608, 24)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 608, 24)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __N.name() : __N,
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_53 = CTD_ANON_53


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_54 (pyxb.binding.basis.complexTypeDefinition):
    """Complex type [anonymous] with content type SIMPLE"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 623, 18)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute N uses Python identifier N
    __N = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'N'), 'N', '__aresysTypes_CTD_ANON_54_N', pyxb.binding.datatypes.int)
    __N._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 626, 24)
    __N._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 626, 24)
    
    N = property(__N.value, __N.set, None, None)

    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_54_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 627, 24)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 627, 24)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __N.name() : __N,
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_54 = CTD_ANON_54


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_55 (pyxb.binding.basis.complexTypeDefinition):
    """Complex type [anonymous] with content type SIMPLE"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 642, 18)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute N uses Python identifier N
    __N = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'N'), 'N', '__aresysTypes_CTD_ANON_55_N', pyxb.binding.datatypes.int)
    __N._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 645, 24)
    __N._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 645, 24)
    
    N = property(__N.value, __N.set, None, None)

    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_55_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 646, 24)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 646, 24)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __N.name() : __N,
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_55 = CTD_ANON_55


# Complex type {aresysTypes}SwathInfoType with content type ELEMENT_ONLY
class SwathInfoType (TreeElementBaseType):
    """Swath general information"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'SwathInfoType')
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 673, 2)
    _ElementMap = TreeElementBaseType._ElementMap.copy()
    _AttributeMap = TreeElementBaseType._AttributeMap.copy()
    # Base type is TreeElementBaseType
    
    # Element Swath uses Python identifier Swath
    __Swath = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'Swath'), 'Swath', '__aresysTypes_SwathInfoType_Swath', False, pyxb.utils.utility.Location('aresysTypes.xsd', 680, 10), )

    
    Swath = property(__Swath.value, __Swath.set, None, 'Swath name')

    
    # Element SwathAcquisitionOrder uses Python identifier SwathAcquisitionOrder
    __SwathAcquisitionOrder = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'SwathAcquisitionOrder'), 'SwathAcquisitionOrder', '__aresysTypes_SwathInfoType_SwathAcquisitionOrder', False, pyxb.utils.utility.Location('aresysTypes.xsd', 692, 10), )

    
    SwathAcquisitionOrder = property(__SwathAcquisitionOrder.value, __SwathAcquisitionOrder.set, None, 'Swath acquisition order')

    
    # Element Polarization uses Python identifier Polarization
    __Polarization = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'Polarization'), 'Polarization', '__aresysTypes_SwathInfoType_Polarization', False, pyxb.utils.utility.Location('aresysTypes.xsd', 704, 10), )

    
    Polarization = property(__Polarization.value, __Polarization.set, None, 'Polarization: H/H, H/V, V/H, V/V')

    
    # Element Rank uses Python identifier Rank
    __Rank = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'Rank'), 'Rank', '__aresysTypes_SwathInfoType_Rank', False, pyxb.utils.utility.Location('aresysTypes.xsd', 709, 10), )

    
    Rank = property(__Rank.value, __Rank.set, None, 'Rank')

    
    # Element RangeDelayBias uses Python identifier RangeDelayBias
    __RangeDelayBias = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'RangeDelayBias'), 'RangeDelayBias', '__aresysTypes_SwathInfoType_RangeDelayBias', False, pyxb.utils.utility.Location('aresysTypes.xsd', 721, 10), )

    
    RangeDelayBias = property(__RangeDelayBias.value, __RangeDelayBias.set, None, 'Range delay bias [s]')

    
    # Element AcquisitionStartTime uses Python identifier AcquisitionStartTime
    __AcquisitionStartTime = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'AcquisitionStartTime'), 'AcquisitionStartTime', '__aresysTypes_SwathInfoType_AcquisitionStartTime', False, pyxb.utils.utility.Location('aresysTypes.xsd', 733, 10), )

    
    AcquisitionStartTime = property(__AcquisitionStartTime.value, __AcquisitionStartTime.set, None, 'Acquisition start time [Utc]')

    
    # Element AzimuthSteeringRateReferenceTime uses Python identifier AzimuthSteeringRateReferenceTime
    __AzimuthSteeringRateReferenceTime = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'AzimuthSteeringRateReferenceTime'), 'AzimuthSteeringRateReferenceTime', '__aresysTypes_SwathInfoType_AzimuthSteeringRateReferenceTime', False, pyxb.utils.utility.Location('aresysTypes.xsd', 745, 10), )

    
    AzimuthSteeringRateReferenceTime = property(__AzimuthSteeringRateReferenceTime.value, __AzimuthSteeringRateReferenceTime.set, None, 'Azimuth antenna steering rate polynomial reference time [s]')

    
    # Element AzimuthSteeringRatePol uses Python identifier AzimuthSteeringRatePol
    __AzimuthSteeringRatePol = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'AzimuthSteeringRatePol'), 'AzimuthSteeringRatePol', '__aresysTypes_SwathInfoType_AzimuthSteeringRatePol', False, pyxb.utils.utility.Location('aresysTypes.xsd', 750, 10), )

    
    AzimuthSteeringRatePol = property(__AzimuthSteeringRatePol.value, __AzimuthSteeringRatePol.set, None, 'Azimuth antenna steering rate polynomial coefficients: const [rad/s], az [rad/s^2], az^2 [rad/s^3]')

    
    # Element AcquisitionPRF uses Python identifier AcquisitionPRF
    __AcquisitionPRF = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'AcquisitionPRF'), 'AcquisitionPRF', '__aresysTypes_SwathInfoType_AcquisitionPRF', False, pyxb.utils.utility.Location('aresysTypes.xsd', 769, 10), )

    
    AcquisitionPRF = property(__AcquisitionPRF.value, __AcquisitionPRF.set, None, 'Acquisition Pulse Repetition Frequency')

    
    # Element EchoesPerBurst uses Python identifier EchoesPerBurst
    __EchoesPerBurst = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'EchoesPerBurst'), 'EchoesPerBurst', '__aresysTypes_SwathInfoType_EchoesPerBurst', False, pyxb.utils.utility.Location('aresysTypes.xsd', 774, 10), )

    
    EchoesPerBurst = property(__EchoesPerBurst.value, __EchoesPerBurst.set, None, 'Number of echoes for each burst')

    
    # Element ChannelDelay uses Python identifier ChannelDelay
    __ChannelDelay = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'ChannelDelay'), 'ChannelDelay', '__aresysTypes_SwathInfoType_ChannelDelay', False, pyxb.utils.utility.Location('aresysTypes.xsd', 779, 10), )

    
    ChannelDelay = property(__ChannelDelay.value, __ChannelDelay.set, None, 'Range channel delay time')

    
    # Element RxGain uses Python identifier RxGain
    __RxGain = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'RxGain'), 'RxGain', '__aresysTypes_SwathInfoType_RxGain', False, pyxb.utils.utility.Location('aresysTypes.xsd', 784, 10), )

    
    RxGain = property(__RxGain.value, __RxGain.set, None, 'Value of the commandable Rx attenuation in the receiver channel')

    
    # Attribute Number inherited from {aresysTypes}TreeElementBaseType
    
    # Attribute Total inherited from {aresysTypes}TreeElementBaseType
    _ElementMap.update({
        __Swath.name() : __Swath,
        __SwathAcquisitionOrder.name() : __SwathAcquisitionOrder,
        __Polarization.name() : __Polarization,
        __Rank.name() : __Rank,
        __RangeDelayBias.name() : __RangeDelayBias,
        __AcquisitionStartTime.name() : __AcquisitionStartTime,
        __AzimuthSteeringRateReferenceTime.name() : __AzimuthSteeringRateReferenceTime,
        __AzimuthSteeringRatePol.name() : __AzimuthSteeringRatePol,
        __AcquisitionPRF.name() : __AcquisitionPRF,
        __EchoesPerBurst.name() : __EchoesPerBurst,
        __ChannelDelay.name() : __ChannelDelay,
        __RxGain.name() : __RxGain
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.SwathInfoType = SwathInfoType
Namespace.addCategoryObject('typeBinding', 'SwathInfoType', SwathInfoType)


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_56 (pyxb.binding.basis.complexTypeDefinition):
    """Swath name"""
    _TypeDefinition = pyxb.binding.datatypes.string
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 684, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.string
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_56_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 687, 18)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 687, 18)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_56 = CTD_ANON_56


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_57 (pyxb.binding.basis.complexTypeDefinition):
    """Swath acquisition order"""
    _TypeDefinition = pyxb.binding.datatypes.unsignedInt
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 696, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.unsignedInt
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_57_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 699, 18)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 699, 18)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_57 = CTD_ANON_57


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_58 (pyxb.binding.basis.complexTypeDefinition):
    """Rank"""
    _TypeDefinition = pyxb.binding.datatypes.int
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 713, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.int
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_58_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 716, 18)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 716, 18)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_58 = CTD_ANON_58


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_59 (pyxb.binding.basis.complexTypeDefinition):
    """Range delay bias [s]"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 725, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_59_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 728, 18)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 728, 18)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_59 = CTD_ANON_59


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_60 (pyxb.binding.basis.complexTypeDefinition):
    """Acquisition start time [Utc]"""
    _TypeDefinition = pyxb.binding.datatypes.string
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 737, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.string
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_60_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 740, 18)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 740, 18)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_60 = CTD_ANON_60


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_61 (pyxb.binding.basis.complexTypeDefinition):
    """Complex type [anonymous] with content type SIMPLE"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 757, 18)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute N uses Python identifier N
    __N = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'N'), 'N', '__aresysTypes_CTD_ANON_61_N', pyxb.binding.datatypes.int)
    __N._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 760, 24)
    __N._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 760, 24)
    
    N = property(__N.value, __N.set, None, None)

    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_61_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 761, 24)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 761, 24)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __N.name() : __N,
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_61 = CTD_ANON_61


# Complex type {aresysTypes}AcquisitionTimelineType with content type ELEMENT_ONLY
class AcquisitionTimelineType (TreeElementBaseType):
    """Acquisition timeline definition"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'AcquisitionTimelineType')
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 793, 2)
    _ElementMap = TreeElementBaseType._ElementMap.copy()
    _AttributeMap = TreeElementBaseType._AttributeMap.copy()
    # Base type is TreeElementBaseType
    
    # Element MissingLines_number uses Python identifier MissingLines_number
    __MissingLines_number = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'MissingLines_number'), 'MissingLines_number', '__aresysTypes_AcquisitionTimelineType_MissingLines_number', False, pyxb.utils.utility.Location('aresysTypes.xsd', 800, 10), )

    
    MissingLines_number = property(__MissingLines_number.value, __MissingLines_number.set, None, 'Number of missing lines')

    
    # Element MissingLines_azimuthtimes uses Python identifier MissingLines_azimuthtimes
    __MissingLines_azimuthtimes = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'MissingLines_azimuthtimes'), 'MissingLines_azimuthtimes', '__aresysTypes_AcquisitionTimelineType_MissingLines_azimuthtimes', False, pyxb.utils.utility.Location('aresysTypes.xsd', 810, 10), )

    
    MissingLines_azimuthtimes = property(__MissingLines_azimuthtimes.value, __MissingLines_azimuthtimes.set, None, 'Azimuth relative times for each missing line')

    
    # Element DuplicatedLines_number uses Python identifier DuplicatedLines_number
    __DuplicatedLines_number = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'DuplicatedLines_number'), 'DuplicatedLines_number', '__aresysTypes_AcquisitionTimelineType_DuplicatedLines_number', False, pyxb.utils.utility.Location('aresysTypes.xsd', 829, 10), )

    
    DuplicatedLines_number = property(__DuplicatedLines_number.value, __DuplicatedLines_number.set, None, 'Number of duplicated lines')

    
    # Element DuplicatedLines_azimuthtimes uses Python identifier DuplicatedLines_azimuthtimes
    __DuplicatedLines_azimuthtimes = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'DuplicatedLines_azimuthtimes'), 'DuplicatedLines_azimuthtimes', '__aresysTypes_AcquisitionTimelineType_DuplicatedLines_azimuthtimes', False, pyxb.utils.utility.Location('aresysTypes.xsd', 839, 10), )

    
    DuplicatedLines_azimuthtimes = property(__DuplicatedLines_azimuthtimes.value, __DuplicatedLines_azimuthtimes.set, None, 'Azimuth relative times for each duplicated line')

    
    # Element PRF_changes_number uses Python identifier PRF_changes_number
    __PRF_changes_number = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'PRF_changes_number'), 'PRF_changes_number', '__aresysTypes_AcquisitionTimelineType_PRF_changes_number', False, pyxb.utils.utility.Location('aresysTypes.xsd', 858, 10), )

    
    PRF_changes_number = property(__PRF_changes_number.value, __PRF_changes_number.set, None, 'Number of PRF changes')

    
    # Element PRF_changes_azimuthtimes uses Python identifier PRF_changes_azimuthtimes
    __PRF_changes_azimuthtimes = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'PRF_changes_azimuthtimes'), 'PRF_changes_azimuthtimes', '__aresysTypes_AcquisitionTimelineType_PRF_changes_azimuthtimes', False, pyxb.utils.utility.Location('aresysTypes.xsd', 868, 10), )

    
    PRF_changes_azimuthtimes = property(__PRF_changes_azimuthtimes.value, __PRF_changes_azimuthtimes.set, None, 'Azimuth relative times for each PRF change')

    
    # Element PRF_changes_values uses Python identifier PRF_changes_values
    __PRF_changes_values = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'PRF_changes_values'), 'PRF_changes_values', '__aresysTypes_AcquisitionTimelineType_PRF_changes_values', False, pyxb.utils.utility.Location('aresysTypes.xsd', 887, 10), )

    
    PRF_changes_values = property(__PRF_changes_values.value, __PRF_changes_values.set, None, 'PRF changes values')

    
    # Element Swst_changes_number uses Python identifier Swst_changes_number
    __Swst_changes_number = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'Swst_changes_number'), 'Swst_changes_number', '__aresysTypes_AcquisitionTimelineType_Swst_changes_number', False, pyxb.utils.utility.Location('aresysTypes.xsd', 906, 10), )

    
    Swst_changes_number = property(__Swst_changes_number.value, __Swst_changes_number.set, None, 'Number of SWST changes')

    
    # Element Swst_changes_azimuthtimes uses Python identifier Swst_changes_azimuthtimes
    __Swst_changes_azimuthtimes = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'Swst_changes_azimuthtimes'), 'Swst_changes_azimuthtimes', '__aresysTypes_AcquisitionTimelineType_Swst_changes_azimuthtimes', False, pyxb.utils.utility.Location('aresysTypes.xsd', 916, 10), )

    
    Swst_changes_azimuthtimes = property(__Swst_changes_azimuthtimes.value, __Swst_changes_azimuthtimes.set, None, 'Azimuth relative times for each SWST change')

    
    # Element Swst_changes_values uses Python identifier Swst_changes_values
    __Swst_changes_values = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'Swst_changes_values'), 'Swst_changes_values', '__aresysTypes_AcquisitionTimelineType_Swst_changes_values', False, pyxb.utils.utility.Location('aresysTypes.xsd', 935, 10), )

    
    Swst_changes_values = property(__Swst_changes_values.value, __Swst_changes_values.set, None, 'SWST changes values')

    
    # Element noise_packets_number uses Python identifier noise_packets_number
    __noise_packets_number = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'noise_packets_number'), 'noise_packets_number', '__aresysTypes_AcquisitionTimelineType_noise_packets_number', False, pyxb.utils.utility.Location('aresysTypes.xsd', 954, 10), )

    
    noise_packets_number = property(__noise_packets_number.value, __noise_packets_number.set, None, 'Number of noise packets')

    
    # Element noise_packets_azimuthtimes uses Python identifier noise_packets_azimuthtimes
    __noise_packets_azimuthtimes = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'noise_packets_azimuthtimes'), 'noise_packets_azimuthtimes', '__aresysTypes_AcquisitionTimelineType_noise_packets_azimuthtimes', False, pyxb.utils.utility.Location('aresysTypes.xsd', 959, 10), )

    
    noise_packets_azimuthtimes = property(__noise_packets_azimuthtimes.value, __noise_packets_azimuthtimes.set, None, 'Azimuth relative times for each noise packet')

    
    # Element Internal_calibration_number uses Python identifier Internal_calibration_number
    __Internal_calibration_number = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'Internal_calibration_number'), 'Internal_calibration_number', '__aresysTypes_AcquisitionTimelineType_Internal_calibration_number', False, pyxb.utils.utility.Location('aresysTypes.xsd', 978, 10), )

    
    Internal_calibration_number = property(__Internal_calibration_number.value, __Internal_calibration_number.set, None, 'Number of internal calibration packets')

    
    # Element Internal_calibration_azimuthtimes uses Python identifier Internal_calibration_azimuthtimes
    __Internal_calibration_azimuthtimes = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'Internal_calibration_azimuthtimes'), 'Internal_calibration_azimuthtimes', '__aresysTypes_AcquisitionTimelineType_Internal_calibration_azimuthtimes', False, pyxb.utils.utility.Location('aresysTypes.xsd', 983, 10), )

    
    Internal_calibration_azimuthtimes = property(__Internal_calibration_azimuthtimes.value, __Internal_calibration_azimuthtimes.set, None, 'Azimuth relative times for each internal calibration packet')

    
    # Element Swl_changes_number uses Python identifier Swl_changes_number
    __Swl_changes_number = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'Swl_changes_number'), 'Swl_changes_number', '__aresysTypes_AcquisitionTimelineType_Swl_changes_number', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1003, 12), )

    
    Swl_changes_number = property(__Swl_changes_number.value, __Swl_changes_number.set, None, 'Number of SWL changes')

    
    # Element Swl_changes_azimuthtimes uses Python identifier Swl_changes_azimuthtimes
    __Swl_changes_azimuthtimes = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'Swl_changes_azimuthtimes'), 'Swl_changes_azimuthtimes', '__aresysTypes_AcquisitionTimelineType_Swl_changes_azimuthtimes', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1008, 12), )

    
    Swl_changes_azimuthtimes = property(__Swl_changes_azimuthtimes.value, __Swl_changes_azimuthtimes.set, None, 'Relative azimuth times for each SWL change')

    
    # Element Swl_changes_values uses Python identifier Swl_changes_values
    __Swl_changes_values = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'Swl_changes_values'), 'Swl_changes_values', '__aresysTypes_AcquisitionTimelineType_Swl_changes_values', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1027, 12), )

    
    Swl_changes_values = property(__Swl_changes_values.value, __Swl_changes_values.set, None, 'SWL changes values')

    
    # Element ChirpPeriod uses Python identifier ChirpPeriod
    __ChirpPeriod = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'ChirpPeriod'), 'ChirpPeriod', '__aresysTypes_AcquisitionTimelineType_ChirpPeriod', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1047, 10), )

    
    ChirpPeriod = property(__ChirpPeriod.value, __ChirpPeriod.set, None, 'Periodic list of chirp indexes related to multi-chirp image acquisition')

    
    # Attribute Number inherited from {aresysTypes}TreeElementBaseType
    
    # Attribute Total inherited from {aresysTypes}TreeElementBaseType
    _ElementMap.update({
        __MissingLines_number.name() : __MissingLines_number,
        __MissingLines_azimuthtimes.name() : __MissingLines_azimuthtimes,
        __DuplicatedLines_number.name() : __DuplicatedLines_number,
        __DuplicatedLines_azimuthtimes.name() : __DuplicatedLines_azimuthtimes,
        __PRF_changes_number.name() : __PRF_changes_number,
        __PRF_changes_azimuthtimes.name() : __PRF_changes_azimuthtimes,
        __PRF_changes_values.name() : __PRF_changes_values,
        __Swst_changes_number.name() : __Swst_changes_number,
        __Swst_changes_azimuthtimes.name() : __Swst_changes_azimuthtimes,
        __Swst_changes_values.name() : __Swst_changes_values,
        __noise_packets_number.name() : __noise_packets_number,
        __noise_packets_azimuthtimes.name() : __noise_packets_azimuthtimes,
        __Internal_calibration_number.name() : __Internal_calibration_number,
        __Internal_calibration_azimuthtimes.name() : __Internal_calibration_azimuthtimes,
        __Swl_changes_number.name() : __Swl_changes_number,
        __Swl_changes_azimuthtimes.name() : __Swl_changes_azimuthtimes,
        __Swl_changes_values.name() : __Swl_changes_values,
        __ChirpPeriod.name() : __ChirpPeriod
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.AcquisitionTimelineType = AcquisitionTimelineType
Namespace.addCategoryObject('typeBinding', 'AcquisitionTimelineType', AcquisitionTimelineType)


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_62 (pyxb.binding.basis.complexTypeDefinition):
    """Complex type [anonymous] with content type SIMPLE"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 817, 18)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_62_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 820, 24)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 820, 24)
    
    unit = property(__unit.value, __unit.set, None, None)

    
    # Attribute N uses Python identifier N
    __N = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'N'), 'N', '__aresysTypes_CTD_ANON_62_N', pyxb.binding.datatypes.unsignedInt)
    __N._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 821, 24)
    __N._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 821, 24)
    
    N = property(__N.value, __N.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit,
        __N.name() : __N
    })
_module_typeBindings.CTD_ANON_62 = CTD_ANON_62


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_63 (pyxb.binding.basis.complexTypeDefinition):
    """Complex type [anonymous] with content type SIMPLE"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 846, 18)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_63_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 849, 24)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 849, 24)
    
    unit = property(__unit.value, __unit.set, None, None)

    
    # Attribute N uses Python identifier N
    __N = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'N'), 'N', '__aresysTypes_CTD_ANON_63_N', pyxb.binding.datatypes.unsignedInt)
    __N._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 850, 24)
    __N._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 850, 24)
    
    N = property(__N.value, __N.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit,
        __N.name() : __N
    })
_module_typeBindings.CTD_ANON_63 = CTD_ANON_63


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_64 (pyxb.binding.basis.complexTypeDefinition):
    """Complex type [anonymous] with content type SIMPLE"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 875, 26)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_64_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 878, 38)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 878, 38)
    
    unit = property(__unit.value, __unit.set, None, None)

    
    # Attribute N uses Python identifier N
    __N = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'N'), 'N', '__aresysTypes_CTD_ANON_64_N', pyxb.binding.datatypes.unsignedInt)
    __N._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 879, 38)
    __N._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 879, 38)
    
    N = property(__N.value, __N.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit,
        __N.name() : __N
    })
_module_typeBindings.CTD_ANON_64 = CTD_ANON_64


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_65 (pyxb.binding.basis.complexTypeDefinition):
    """Complex type [anonymous] with content type SIMPLE"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 894, 26)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_65_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 897, 38)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 897, 38)
    
    unit = property(__unit.value, __unit.set, None, None)

    
    # Attribute N uses Python identifier N
    __N = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'N'), 'N', '__aresysTypes_CTD_ANON_65_N', pyxb.binding.datatypes.unsignedInt)
    __N._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 898, 38)
    __N._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 898, 38)
    
    N = property(__N.value, __N.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit,
        __N.name() : __N
    })
_module_typeBindings.CTD_ANON_65 = CTD_ANON_65


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_66 (pyxb.binding.basis.complexTypeDefinition):
    """Complex type [anonymous] with content type SIMPLE"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 923, 18)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_66_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 926, 24)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 926, 24)
    
    unit = property(__unit.value, __unit.set, None, None)

    
    # Attribute N uses Python identifier N
    __N = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'N'), 'N', '__aresysTypes_CTD_ANON_66_N', pyxb.binding.datatypes.unsignedInt)
    __N._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 927, 24)
    __N._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 927, 24)
    
    N = property(__N.value, __N.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit,
        __N.name() : __N
    })
_module_typeBindings.CTD_ANON_66 = CTD_ANON_66


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_67 (pyxb.binding.basis.complexTypeDefinition):
    """Complex type [anonymous] with content type SIMPLE"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 942, 18)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_67_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 945, 24)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 945, 24)
    
    unit = property(__unit.value, __unit.set, None, None)

    
    # Attribute N uses Python identifier N
    __N = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'N'), 'N', '__aresysTypes_CTD_ANON_67_N', pyxb.binding.datatypes.unsignedInt)
    __N._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 946, 24)
    __N._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 946, 24)
    
    N = property(__N.value, __N.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit,
        __N.name() : __N
    })
_module_typeBindings.CTD_ANON_67 = CTD_ANON_67


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_68 (pyxb.binding.basis.complexTypeDefinition):
    """Complex type [anonymous] with content type SIMPLE"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 966, 18)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute N uses Python identifier N
    __N = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'N'), 'N', '__aresysTypes_CTD_ANON_68_N', pyxb.binding.datatypes.unsignedInt)
    __N._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 969, 24)
    __N._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 969, 24)
    
    N = property(__N.value, __N.set, None, None)

    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_68_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 970, 24)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 970, 24)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __N.name() : __N,
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_68 = CTD_ANON_68


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_69 (pyxb.binding.basis.complexTypeDefinition):
    """Complex type [anonymous] with content type SIMPLE"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 990, 18)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute N uses Python identifier N
    __N = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'N'), 'N', '__aresysTypes_CTD_ANON_69_N', pyxb.binding.datatypes.unsignedInt)
    __N._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 993, 24)
    __N._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 993, 24)
    
    N = property(__N.value, __N.set, None, None)

    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_69_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 994, 24)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 994, 24)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __N.name() : __N,
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_69 = CTD_ANON_69


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_70 (pyxb.binding.basis.complexTypeDefinition):
    """Complex type [anonymous] with content type SIMPLE"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1015, 20)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_70_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1018, 26)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1018, 26)
    
    unit = property(__unit.value, __unit.set, None, None)

    
    # Attribute N uses Python identifier N
    __N = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'N'), 'N', '__aresysTypes_CTD_ANON_70_N', pyxb.binding.datatypes.unsignedInt)
    __N._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1019, 26)
    __N._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1019, 26)
    
    N = property(__N.value, __N.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit,
        __N.name() : __N
    })
_module_typeBindings.CTD_ANON_70 = CTD_ANON_70


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_71 (pyxb.binding.basis.complexTypeDefinition):
    """Complex type [anonymous] with content type SIMPLE"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1034, 20)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_71_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1037, 26)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1037, 26)
    
    unit = property(__unit.value, __unit.set, None, None)

    
    # Attribute N uses Python identifier N
    __N = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'N'), 'N', '__aresysTypes_CTD_ANON_71_N', pyxb.binding.datatypes.unsignedInt)
    __N._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1038, 26)
    __N._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1038, 26)
    
    N = property(__N.value, __N.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit,
        __N.name() : __N
    })
_module_typeBindings.CTD_ANON_71 = CTD_ANON_71


# Complex type {aresysTypes}RasterInfoType with content type ELEMENT_ONLY
class RasterInfoType (TreeElementBaseType):
    """Information regarding binary file format and time coordinates of the image"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'RasterInfoType')
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1056, 2)
    _ElementMap = TreeElementBaseType._ElementMap.copy()
    _AttributeMap = TreeElementBaseType._AttributeMap.copy()
    # Base type is TreeElementBaseType
    
    # Element FileName uses Python identifier FileName
    __FileName = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'FileName'), 'FileName', '__aresysTypes_RasterInfoType_FileName', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1063, 10), )

    
    FileName = property(__FileName.value, __FileName.set, None, 'Name of the associated binary file')

    
    # Element Lines uses Python identifier Lines
    __Lines = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'Lines'), 'Lines', '__aresysTypes_RasterInfoType_Lines', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1068, 10), )

    
    Lines = property(__Lines.value, __Lines.set, None, 'Total number of lines (azimuth) of the image')

    
    # Element Samples uses Python identifier Samples
    __Samples = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'Samples'), 'Samples', '__aresysTypes_RasterInfoType_Samples', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1073, 10), )

    
    Samples = property(__Samples.value, __Samples.set, None, 'Total number of samples (range) of the image')

    
    # Element HeaderOffsetBytes uses Python identifier HeaderOffsetBytes
    __HeaderOffsetBytes = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'HeaderOffsetBytes'), 'HeaderOffsetBytes', '__aresysTypes_RasterInfoType_HeaderOffsetBytes', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1078, 10), )

    
    HeaderOffsetBytes = property(__HeaderOffsetBytes.value, __HeaderOffsetBytes.set, None, 'Number of bytes at the beginning of the file containing the header information')

    
    # Element RowPrefixBytes uses Python identifier RowPrefixBytes
    __RowPrefixBytes = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'RowPrefixBytes'), 'RowPrefixBytes', '__aresysTypes_RasterInfoType_RowPrefixBytes', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1083, 10), )

    
    RowPrefixBytes = property(__RowPrefixBytes.value, __RowPrefixBytes.set, None, 'Number of bytes at the beginning of each line containing header information')

    
    # Element ByteOrder uses Python identifier ByteOrder
    __ByteOrder = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'ByteOrder'), 'ByteOrder', '__aresysTypes_RasterInfoType_ByteOrder', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1088, 10), )

    
    ByteOrder = property(__ByteOrder.value, __ByteOrder.set, None, 'Endianity: BIGENDIAN or LITTLEENDIAN')

    
    # Element CellType uses Python identifier CellType
    __CellType = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'CellType'), 'CellType', '__aresysTypes_RasterInfoType_CellType', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1093, 10), )

    
    CellType = property(__CellType.value, __CellType.set, None, 'Byte format type: FLOAT_COMPLEX, FLOAT32, DOUBLE_COMPLEX, FLOAT64, INT16, SHORT_COMPLEX, INT32, INT_COMPLEX, INT8, INT8_COMPLEX')

    
    # Element LinesStep uses Python identifier LinesStep
    __LinesStep = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'LinesStep'), 'LinesStep', '__aresysTypes_RasterInfoType_LinesStep', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1098, 10), )

    
    LinesStep = property(__LinesStep.value, __LinesStep.set, None, 'Azimuth sampling step [s]')

    
    # Element SamplesStep uses Python identifier SamplesStep
    __SamplesStep = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'SamplesStep'), 'SamplesStep', '__aresysTypes_RasterInfoType_SamplesStep', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1103, 10), )

    
    SamplesStep = property(__SamplesStep.value, __SamplesStep.set, None, 'Range sampling step [s]')

    
    # Element LinesStart uses Python identifier LinesStart
    __LinesStart = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'LinesStart'), 'LinesStart', '__aresysTypes_RasterInfoType_LinesStart', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1108, 10), )

    
    LinesStart = property(__LinesStart.value, __LinesStart.set, None, 'Azimuth absolute start time [Utc]')

    
    # Element SamplesStart uses Python identifier SamplesStart
    __SamplesStart = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'SamplesStart'), 'SamplesStart', '__aresysTypes_RasterInfoType_SamplesStart', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1113, 10), )

    
    SamplesStart = property(__SamplesStart.value, __SamplesStart.set, None, 'Range absolute start time [s]')

    
    # Element RasterFormat uses Python identifier RasterFormat
    __RasterFormat = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'RasterFormat'), 'RasterFormat', '__aresysTypes_RasterInfoType_RasterFormat', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1118, 10), )

    
    RasterFormat = property(__RasterFormat.value, __RasterFormat.set, None, 'Raster Format of the Data ( Default value is ARESYS_RASTER )')

    
    # Element InvalidValue uses Python identifier InvalidValue
    __InvalidValue = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'InvalidValue'), 'InvalidValue', '__aresysTypes_RasterInfoType_InvalidValue', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1123, 10), )

    
    InvalidValue = property(__InvalidValue.value, __InvalidValue.set, None, None)

    
    # Attribute Number inherited from {aresysTypes}TreeElementBaseType
    
    # Attribute Total inherited from {aresysTypes}TreeElementBaseType
    _ElementMap.update({
        __FileName.name() : __FileName,
        __Lines.name() : __Lines,
        __Samples.name() : __Samples,
        __HeaderOffsetBytes.name() : __HeaderOffsetBytes,
        __RowPrefixBytes.name() : __RowPrefixBytes,
        __ByteOrder.name() : __ByteOrder,
        __CellType.name() : __CellType,
        __LinesStep.name() : __LinesStep,
        __SamplesStep.name() : __SamplesStep,
        __LinesStart.name() : __LinesStart,
        __SamplesStart.name() : __SamplesStart,
        __RasterFormat.name() : __RasterFormat,
        __InvalidValue.name() : __InvalidValue
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.RasterInfoType = RasterInfoType
Namespace.addCategoryObject('typeBinding', 'RasterInfoType', RasterInfoType)


# Complex type {aresysTypes}SamplingConstantsType with content type ELEMENT_ONLY
class SamplingConstantsType (TreeElementBaseType):
    """Bandwidths and sampling frequencies"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'SamplingConstantsType')
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1128, 2)
    _ElementMap = TreeElementBaseType._ElementMap.copy()
    _AttributeMap = TreeElementBaseType._AttributeMap.copy()
    # Base type is TreeElementBaseType
    
    # Element frg_hz uses Python identifier frg_hz
    __frg_hz = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'frg_hz'), 'frg_hz', '__aresysTypes_SamplingConstantsType_frg_hz', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1135, 10), )

    
    frg_hz = property(__frg_hz.value, __frg_hz.set, None, 'Range sampling frequency [Hz]')

    
    # Element Brg_hz uses Python identifier Brg_hz
    __Brg_hz = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'Brg_hz'), 'Brg_hz', '__aresysTypes_SamplingConstantsType_Brg_hz', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1147, 10), )

    
    Brg_hz = property(__Brg_hz.value, __Brg_hz.set, None, 'Range bandwidth [Hz]')

    
    # Element faz_hz uses Python identifier faz_hz
    __faz_hz = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'faz_hz'), 'faz_hz', '__aresysTypes_SamplingConstantsType_faz_hz', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1159, 10), )

    
    faz_hz = property(__faz_hz.value, __faz_hz.set, None, 'Azimuth sampling frequency [Hz]')

    
    # Element Baz_hz uses Python identifier Baz_hz
    __Baz_hz = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'Baz_hz'), 'Baz_hz', '__aresysTypes_SamplingConstantsType_Baz_hz', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1171, 10), )

    
    Baz_hz = property(__Baz_hz.value, __Baz_hz.set, None, 'Azimuth bandwidth [Hz]')

    
    # Attribute Number inherited from {aresysTypes}TreeElementBaseType
    
    # Attribute Total inherited from {aresysTypes}TreeElementBaseType
    _ElementMap.update({
        __frg_hz.name() : __frg_hz,
        __Brg_hz.name() : __Brg_hz,
        __faz_hz.name() : __faz_hz,
        __Baz_hz.name() : __Baz_hz
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.SamplingConstantsType = SamplingConstantsType
Namespace.addCategoryObject('typeBinding', 'SamplingConstantsType', SamplingConstantsType)


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_72 (pyxb.binding.basis.complexTypeDefinition):
    """Range sampling frequency [Hz]"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1139, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_72_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1142, 18)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1142, 18)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_72 = CTD_ANON_72


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_73 (pyxb.binding.basis.complexTypeDefinition):
    """Range bandwidth [Hz]"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1151, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_73_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1154, 18)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1154, 18)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_73 = CTD_ANON_73


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_74 (pyxb.binding.basis.complexTypeDefinition):
    """Azimuth sampling frequency [Hz]"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1163, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_74_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1166, 18)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1166, 18)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_74 = CTD_ANON_74


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_75 (pyxb.binding.basis.complexTypeDefinition):
    """Azimuth bandwidth [Hz]"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1175, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_75_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1178, 18)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1178, 18)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_75 = CTD_ANON_75


# Complex type {aresysTypes}DataStatisticsType with content type ELEMENT_ONLY
class DataStatisticsType (TreeElementBaseType):
    """Statistics computed from data"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'DataStatisticsType')
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1187, 2)
    _ElementMap = TreeElementBaseType._ElementMap.copy()
    _AttributeMap = TreeElementBaseType._AttributeMap.copy()
    # Base type is TreeElementBaseType
    
    # Element NumSamples uses Python identifier NumSamples
    __NumSamples = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'NumSamples'), 'NumSamples', '__aresysTypes_DataStatisticsType_NumSamples', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1194, 10), )

    
    NumSamples = property(__NumSamples.value, __NumSamples.set, None, 'Number of samples analyzed')

    
    # Element MaxI uses Python identifier MaxI
    __MaxI = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'MaxI'), 'MaxI', '__aresysTypes_DataStatisticsType_MaxI', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1206, 10), )

    
    MaxI = property(__MaxI.value, __MaxI.set, None, 'Max of I (real) samples')

    
    # Element MinI uses Python identifier MinI
    __MinI = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'MinI'), 'MinI', '__aresysTypes_DataStatisticsType_MinI', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1218, 10), )

    
    MinI = property(__MinI.value, __MinI.set, None, 'Min of I (real) samples')

    
    # Element MaxQ uses Python identifier MaxQ
    __MaxQ = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'MaxQ'), 'MaxQ', '__aresysTypes_DataStatisticsType_MaxQ', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1230, 10), )

    
    MaxQ = property(__MaxQ.value, __MaxQ.set, None, 'Max of Q (imaginary) samples')

    
    # Element MinQ uses Python identifier MinQ
    __MinQ = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'MinQ'), 'MinQ', '__aresysTypes_DataStatisticsType_MinQ', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1242, 10), )

    
    MinQ = property(__MinQ.value, __MinQ.set, None, 'Min of Q (imaginary) samples')

    
    # Element SumI uses Python identifier SumI
    __SumI = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'SumI'), 'SumI', '__aresysTypes_DataStatisticsType_SumI', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1254, 10), )

    
    SumI = property(__SumI.value, __SumI.set, None, 'Sum of I (real) samples')

    
    # Element SumQ uses Python identifier SumQ
    __SumQ = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'SumQ'), 'SumQ', '__aresysTypes_DataStatisticsType_SumQ', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1266, 10), )

    
    SumQ = property(__SumQ.value, __SumQ.set, None, 'Sum of Q (imaginary) samples')

    
    # Element Sum2I uses Python identifier Sum2I
    __Sum2I = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'Sum2I'), 'Sum2I', '__aresysTypes_DataStatisticsType_Sum2I', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1278, 10), )

    
    Sum2I = property(__Sum2I.value, __Sum2I.set, None, 'Square Sum of I (real) samples')

    
    # Element Sum2Q uses Python identifier Sum2Q
    __Sum2Q = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'Sum2Q'), 'Sum2Q', '__aresysTypes_DataStatisticsType_Sum2Q', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1290, 10), )

    
    Sum2Q = property(__Sum2Q.value, __Sum2Q.set, None, 'Square Sum of Q (imaginary) samples')

    
    # Element StdDevI uses Python identifier StdDevI
    __StdDevI = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'StdDevI'), 'StdDevI', '__aresysTypes_DataStatisticsType_StdDevI', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1302, 10), )

    
    StdDevI = property(__StdDevI.value, __StdDevI.set, None, 'Standard Deviation of I (real) samples')

    
    # Element StdDevQ uses Python identifier StdDevQ
    __StdDevQ = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'StdDevQ'), 'StdDevQ', '__aresysTypes_DataStatisticsType_StdDevQ', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1314, 10), )

    
    StdDevQ = property(__StdDevQ.value, __StdDevQ.set, None, 'Standard Deviation of Q (imaginary) samples')

    
    # Element StatisticsList uses Python identifier StatisticsList
    __StatisticsList = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'StatisticsList'), 'StatisticsList', '__aresysTypes_DataStatisticsType_StatisticsList', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1326, 10), )

    
    StatisticsList = property(__StatisticsList.value, __StatisticsList.set, None, None)

    
    # Attribute Number inherited from {aresysTypes}TreeElementBaseType
    
    # Attribute Total inherited from {aresysTypes}TreeElementBaseType
    _ElementMap.update({
        __NumSamples.name() : __NumSamples,
        __MaxI.name() : __MaxI,
        __MinI.name() : __MinI,
        __MaxQ.name() : __MaxQ,
        __MinQ.name() : __MinQ,
        __SumI.name() : __SumI,
        __SumQ.name() : __SumQ,
        __Sum2I.name() : __Sum2I,
        __Sum2Q.name() : __Sum2Q,
        __StdDevI.name() : __StdDevI,
        __StdDevQ.name() : __StdDevQ,
        __StatisticsList.name() : __StatisticsList
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.DataStatisticsType = DataStatisticsType
Namespace.addCategoryObject('typeBinding', 'DataStatisticsType', DataStatisticsType)


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_76 (pyxb.binding.basis.complexTypeDefinition):
    """Number of samples analyzed"""
    _TypeDefinition = pyxb.binding.datatypes.unsignedLong
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1198, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.unsignedLong
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_76_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1201, 18)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1201, 18)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_76 = CTD_ANON_76


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_77 (pyxb.binding.basis.complexTypeDefinition):
    """Max of I (real) samples"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1210, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_77_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1213, 18)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1213, 18)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_77 = CTD_ANON_77


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_78 (pyxb.binding.basis.complexTypeDefinition):
    """Min of I (real) samples"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1222, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_78_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1225, 18)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1225, 18)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_78 = CTD_ANON_78


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_79 (pyxb.binding.basis.complexTypeDefinition):
    """Max of Q (imaginary) samples"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1234, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_79_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1237, 18)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1237, 18)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_79 = CTD_ANON_79


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_80 (pyxb.binding.basis.complexTypeDefinition):
    """Min of Q (imaginary) samples"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1246, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_80_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1249, 18)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1249, 18)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_80 = CTD_ANON_80


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_81 (pyxb.binding.basis.complexTypeDefinition):
    """Sum of I (real) samples"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1258, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_81_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1261, 18)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1261, 18)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_81 = CTD_ANON_81


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_82 (pyxb.binding.basis.complexTypeDefinition):
    """Sum of Q (imaginary) samples"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1270, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_82_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1273, 18)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1273, 18)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_82 = CTD_ANON_82


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_83 (pyxb.binding.basis.complexTypeDefinition):
    """Square Sum of I (real) samples"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1282, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_83_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1285, 18)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1285, 18)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_83 = CTD_ANON_83


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_84 (pyxb.binding.basis.complexTypeDefinition):
    """Square Sum of Q (imaginary) samples"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1294, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_84_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1297, 18)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1297, 18)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_84 = CTD_ANON_84


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_85 (pyxb.binding.basis.complexTypeDefinition):
    """Standard Deviation of I (real) samples"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1306, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_85_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1309, 18)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1309, 18)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_85 = CTD_ANON_85


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_86 (pyxb.binding.basis.complexTypeDefinition):
    """Standard Deviation of Q (imaginary) samples"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1318, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_86_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1321, 18)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1321, 18)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_86 = CTD_ANON_86


# Complex type {aresysTypes}DataBlockStatisticsType with content type ELEMENT_ONLY
class DataBlockStatisticsType (TreeElementBaseType):
    """Statistics computed from data"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'DataBlockStatisticsType')
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1337, 2)
    _ElementMap = TreeElementBaseType._ElementMap.copy()
    _AttributeMap = TreeElementBaseType._AttributeMap.copy()
    # Base type is TreeElementBaseType
    
    # Element NumSamples uses Python identifier NumSamples
    __NumSamples = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'NumSamples'), 'NumSamples', '__aresysTypes_DataBlockStatisticsType_NumSamples', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1344, 10), )

    
    NumSamples = property(__NumSamples.value, __NumSamples.set, None, 'Number of samples analyzed')

    
    # Element MaxI uses Python identifier MaxI
    __MaxI = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'MaxI'), 'MaxI', '__aresysTypes_DataBlockStatisticsType_MaxI', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1356, 10), )

    
    MaxI = property(__MaxI.value, __MaxI.set, None, 'Max of I (real) samples')

    
    # Element MinI uses Python identifier MinI
    __MinI = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'MinI'), 'MinI', '__aresysTypes_DataBlockStatisticsType_MinI', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1368, 10), )

    
    MinI = property(__MinI.value, __MinI.set, None, 'Min of I (real) samples')

    
    # Element MaxQ uses Python identifier MaxQ
    __MaxQ = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'MaxQ'), 'MaxQ', '__aresysTypes_DataBlockStatisticsType_MaxQ', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1380, 10), )

    
    MaxQ = property(__MaxQ.value, __MaxQ.set, None, 'Max of Q (imaginary) samples')

    
    # Element MinQ uses Python identifier MinQ
    __MinQ = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'MinQ'), 'MinQ', '__aresysTypes_DataBlockStatisticsType_MinQ', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1392, 10), )

    
    MinQ = property(__MinQ.value, __MinQ.set, None, 'Min of Q (imaginary) samples')

    
    # Element SumI uses Python identifier SumI
    __SumI = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'SumI'), 'SumI', '__aresysTypes_DataBlockStatisticsType_SumI', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1404, 10), )

    
    SumI = property(__SumI.value, __SumI.set, None, 'Sum of I (real) samples')

    
    # Element SumQ uses Python identifier SumQ
    __SumQ = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'SumQ'), 'SumQ', '__aresysTypes_DataBlockStatisticsType_SumQ', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1416, 10), )

    
    SumQ = property(__SumQ.value, __SumQ.set, None, 'Sum of Q (imaginary) samples')

    
    # Element Sum2I uses Python identifier Sum2I
    __Sum2I = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'Sum2I'), 'Sum2I', '__aresysTypes_DataBlockStatisticsType_Sum2I', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1428, 10), )

    
    Sum2I = property(__Sum2I.value, __Sum2I.set, None, 'Square Sum of I (real) samples')

    
    # Element Sum2Q uses Python identifier Sum2Q
    __Sum2Q = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'Sum2Q'), 'Sum2Q', '__aresysTypes_DataBlockStatisticsType_Sum2Q', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1440, 10), )

    
    Sum2Q = property(__Sum2Q.value, __Sum2Q.set, None, 'Square Sum of Q (imaginary) samples')

    
    # Attribute Number inherited from {aresysTypes}TreeElementBaseType
    
    # Attribute Total inherited from {aresysTypes}TreeElementBaseType
    
    # Attribute lineStart uses Python identifier lineStart
    __lineStart = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'lineStart'), 'lineStart', '__aresysTypes_DataBlockStatisticsType_lineStart', pyxb.binding.datatypes.unsignedLong, required=True)
    __lineStart._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1453, 8)
    __lineStart._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1453, 8)
    
    lineStart = property(__lineStart.value, __lineStart.set, None, None)

    
    # Attribute lineStop uses Python identifier lineStop
    __lineStop = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'lineStop'), 'lineStop', '__aresysTypes_DataBlockStatisticsType_lineStop', pyxb.binding.datatypes.unsignedLong, required=True)
    __lineStop._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1454, 8)
    __lineStop._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1454, 8)
    
    lineStop = property(__lineStop.value, __lineStop.set, None, None)

    _ElementMap.update({
        __NumSamples.name() : __NumSamples,
        __MaxI.name() : __MaxI,
        __MinI.name() : __MinI,
        __MaxQ.name() : __MaxQ,
        __MinQ.name() : __MinQ,
        __SumI.name() : __SumI,
        __SumQ.name() : __SumQ,
        __Sum2I.name() : __Sum2I,
        __Sum2Q.name() : __Sum2Q
    })
    _AttributeMap.update({
        __lineStart.name() : __lineStart,
        __lineStop.name() : __lineStop
    })
_module_typeBindings.DataBlockStatisticsType = DataBlockStatisticsType
Namespace.addCategoryObject('typeBinding', 'DataBlockStatisticsType', DataBlockStatisticsType)


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_87 (pyxb.binding.basis.complexTypeDefinition):
    """Number of samples analyzed"""
    _TypeDefinition = pyxb.binding.datatypes.unsignedLong
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1348, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.unsignedLong
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_87_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1351, 18)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1351, 18)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_87 = CTD_ANON_87


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_88 (pyxb.binding.basis.complexTypeDefinition):
    """Max of I (real) samples"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1360, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_88_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1363, 18)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1363, 18)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_88 = CTD_ANON_88


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_89 (pyxb.binding.basis.complexTypeDefinition):
    """Min of I (real) samples"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1372, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_89_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1375, 18)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1375, 18)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_89 = CTD_ANON_89


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_90 (pyxb.binding.basis.complexTypeDefinition):
    """Max of Q (imaginary) samples"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1384, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_90_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1387, 18)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1387, 18)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_90 = CTD_ANON_90


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_91 (pyxb.binding.basis.complexTypeDefinition):
    """Min of Q (imaginary) samples"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1396, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_91_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1399, 18)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1399, 18)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_91 = CTD_ANON_91


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_92 (pyxb.binding.basis.complexTypeDefinition):
    """Sum of I (real) samples"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1408, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_92_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1411, 18)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1411, 18)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_92 = CTD_ANON_92


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_93 (pyxb.binding.basis.complexTypeDefinition):
    """Sum of Q (imaginary) samples"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1420, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_93_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1423, 18)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1423, 18)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_93 = CTD_ANON_93


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_94 (pyxb.binding.basis.complexTypeDefinition):
    """Square Sum of I (real) samples"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1432, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_94_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1435, 18)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1435, 18)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_94 = CTD_ANON_94


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_95 (pyxb.binding.basis.complexTypeDefinition):
    """Square Sum of Q (imaginary) samples"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1444, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_95_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1447, 18)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1447, 18)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_95 = CTD_ANON_95


# Complex type {aresysTypes}polyType with content type ELEMENT_ONLY
class polyType (TreeElementBaseType):
    """Polynomial parametrization of the geometry parameter"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'polyType')
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1458, 2)
    _ElementMap = TreeElementBaseType._ElementMap.copy()
    _AttributeMap = TreeElementBaseType._AttributeMap.copy()
    # Base type is TreeElementBaseType
    
    # Element pol uses Python identifier pol
    __pol = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'pol'), 'pol', '__aresysTypes_polyType_pol', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1465, 10), )

    
    pol = property(__pol.value, __pol.set, None, 'Polynomial coefficients: const, rg, az, az*rg, rg^2, rg^3, rg^4 [Optional: rg^5 rg^6 .... rg^N]')

    
    # Element trg0_s uses Python identifier trg0_s
    __trg0_s = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'trg0_s'), 'trg0_s', '__aresysTypes_polyType_trg0_s', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1484, 10), )

    
    trg0_s = property(__trg0_s.value, __trg0_s.set, None, 'Polynomial range reference time [s]')

    
    # Element taz0_Utc uses Python identifier taz0_Utc
    __taz0_Utc = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'taz0_Utc'), 'taz0_Utc', '__aresysTypes_polyType_taz0_Utc', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1496, 10), )

    
    taz0_Utc = property(__taz0_Utc.value, __taz0_Utc.set, None, 'Polynomial azimuth reference time [Utc]')

    
    # Attribute Number inherited from {aresysTypes}TreeElementBaseType
    
    # Attribute Total inherited from {aresysTypes}TreeElementBaseType
    _ElementMap.update({
        __pol.name() : __pol,
        __trg0_s.name() : __trg0_s,
        __taz0_Utc.name() : __taz0_Utc
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.polyType = polyType
Namespace.addCategoryObject('typeBinding', 'polyType', polyType)


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_96 (pyxb.binding.basis.complexTypeDefinition):
    """Complex type [anonymous] with content type SIMPLE"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1472, 18)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute N uses Python identifier N
    __N = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'N'), 'N', '__aresysTypes_CTD_ANON_96_N', pyxb.binding.datatypes.int)
    __N._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1475, 24)
    __N._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1475, 24)
    
    N = property(__N.value, __N.set, None, None)

    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_96_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1476, 24)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1476, 24)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __N.name() : __N,
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_96 = CTD_ANON_96


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_97 (pyxb.binding.basis.complexTypeDefinition):
    """Polynomial range reference time [s]"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1488, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_97_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1491, 18)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1491, 18)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_97 = CTD_ANON_97


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_98 (pyxb.binding.basis.complexTypeDefinition):
    """Polynomial azimuth reference time [Utc]"""
    _TypeDefinition = pyxb.binding.datatypes.string
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1500, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.string
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_98_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1503, 18)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1503, 18)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_98 = CTD_ANON_98


# Complex type {aresysTypes}polyCoregType with content type ELEMENT_ONLY
class polyCoregType (TreeElementBaseType):
    """Polynomial parametrization of the coregistration parameters"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'polyCoregType')
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1512, 2)
    _ElementMap = TreeElementBaseType._ElementMap.copy()
    _AttributeMap = TreeElementBaseType._AttributeMap.copy()
    # Base type is TreeElementBaseType
    
    # Element polRg uses Python identifier polRg
    __polRg = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'polRg'), 'polRg', '__aresysTypes_polyCoregType_polRg', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1519, 10), )

    
    polRg = property(__polRg.value, __polRg.set, None, 'Polynomial coefficients of Range')

    
    # Element polAz uses Python identifier polAz
    __polAz = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'polAz'), 'polAz', '__aresysTypes_polyCoregType_polAz', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1538, 10), )

    
    polAz = property(__polAz.value, __polAz.set, None, 'Polynomial coefficients of Azimuth')

    
    # Element trg0_s uses Python identifier trg0_s
    __trg0_s = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'trg0_s'), 'trg0_s', '__aresysTypes_polyCoregType_trg0_s', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1557, 10), )

    
    trg0_s = property(__trg0_s.value, __trg0_s.set, None, 'Polynomial range reference time [s]')

    
    # Element taz0_Utc uses Python identifier taz0_Utc
    __taz0_Utc = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'taz0_Utc'), 'taz0_Utc', '__aresysTypes_polyCoregType_taz0_Utc', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1569, 10), )

    
    taz0_Utc = property(__taz0_Utc.value, __taz0_Utc.set, None, 'Polynomial azimuth reference time [Utc]')

    
    # Attribute Number inherited from {aresysTypes}TreeElementBaseType
    
    # Attribute Total inherited from {aresysTypes}TreeElementBaseType
    _ElementMap.update({
        __polRg.name() : __polRg,
        __polAz.name() : __polAz,
        __trg0_s.name() : __trg0_s,
        __taz0_Utc.name() : __taz0_Utc
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.polyCoregType = polyCoregType
Namespace.addCategoryObject('typeBinding', 'polyCoregType', polyCoregType)


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_99 (pyxb.binding.basis.complexTypeDefinition):
    """Complex type [anonymous] with content type SIMPLE"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1526, 18)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute N uses Python identifier N
    __N = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'N'), 'N', '__aresysTypes_CTD_ANON_99_N', pyxb.binding.datatypes.int)
    __N._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1529, 24)
    __N._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1529, 24)
    
    N = property(__N.value, __N.set, None, None)

    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_99_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1530, 24)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1530, 24)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __N.name() : __N,
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_99 = CTD_ANON_99


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_100 (pyxb.binding.basis.complexTypeDefinition):
    """Complex type [anonymous] with content type SIMPLE"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1545, 18)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute N uses Python identifier N
    __N = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'N'), 'N', '__aresysTypes_CTD_ANON_100_N', pyxb.binding.datatypes.int)
    __N._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1548, 24)
    __N._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1548, 24)
    
    N = property(__N.value, __N.set, None, None)

    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_100_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1549, 24)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1549, 24)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __N.name() : __N,
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_100 = CTD_ANON_100


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_101 (pyxb.binding.basis.complexTypeDefinition):
    """Polynomial range reference time [s]"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1561, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_101_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1564, 18)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1564, 18)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_101 = CTD_ANON_101


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_102 (pyxb.binding.basis.complexTypeDefinition):
    """Polynomial azimuth reference time [Utc]"""
    _TypeDefinition = pyxb.binding.datatypes.string
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1573, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.string
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_102_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1576, 18)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1576, 18)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_102 = CTD_ANON_102


# Complex type {aresysTypes}GroundCornersPointsType with content type ELEMENT_ONLY
class GroundCornersPointsType (TreeElementBaseType):
    """Complex type {aresysTypes}GroundCornersPointsType with content type ELEMENT_ONLY"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'GroundCornersPointsType')
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1585, 2)
    _ElementMap = TreeElementBaseType._ElementMap.copy()
    _AttributeMap = TreeElementBaseType._AttributeMap.copy()
    # Base type is TreeElementBaseType
    
    # Element EastingGridSize uses Python identifier EastingGridSize
    __EastingGridSize = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'EastingGridSize'), 'EastingGridSize', '__aresysTypes_GroundCornersPointsType_EastingGridSize', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1589, 10), )

    
    EastingGridSize = property(__EastingGridSize.value, __EastingGridSize.set, None, None)

    
    # Element NorthingGridSize uses Python identifier NorthingGridSize
    __NorthingGridSize = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'NorthingGridSize'), 'NorthingGridSize', '__aresysTypes_GroundCornersPointsType_NorthingGridSize', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1598, 10), )

    
    NorthingGridSize = property(__NorthingGridSize.value, __NorthingGridSize.set, None, None)

    
    # Element NorthWest uses Python identifier NorthWest
    __NorthWest = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'NorthWest'), 'NorthWest', '__aresysTypes_GroundCornersPointsType_NorthWest', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1607, 10), )

    
    NorthWest = property(__NorthWest.value, __NorthWest.set, None, None)

    
    # Element NorthEast uses Python identifier NorthEast
    __NorthEast = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'NorthEast'), 'NorthEast', '__aresysTypes_GroundCornersPointsType_NorthEast', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1614, 10), )

    
    NorthEast = property(__NorthEast.value, __NorthEast.set, None, None)

    
    # Element SouthWest uses Python identifier SouthWest
    __SouthWest = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'SouthWest'), 'SouthWest', '__aresysTypes_GroundCornersPointsType_SouthWest', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1621, 10), )

    
    SouthWest = property(__SouthWest.value, __SouthWest.set, None, None)

    
    # Element SouthEast uses Python identifier SouthEast
    __SouthEast = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'SouthEast'), 'SouthEast', '__aresysTypes_GroundCornersPointsType_SouthEast', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1628, 10), )

    
    SouthEast = property(__SouthEast.value, __SouthEast.set, None, None)

    
    # Element Center uses Python identifier Center
    __Center = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'Center'), 'Center', '__aresysTypes_GroundCornersPointsType_Center', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1635, 10), )

    
    Center = property(__Center.value, __Center.set, None, None)

    
    # Attribute Number inherited from {aresysTypes}TreeElementBaseType
    
    # Attribute Total inherited from {aresysTypes}TreeElementBaseType
    _ElementMap.update({
        __EastingGridSize.name() : __EastingGridSize,
        __NorthingGridSize.name() : __NorthingGridSize,
        __NorthWest.name() : __NorthWest,
        __NorthEast.name() : __NorthEast,
        __SouthWest.name() : __SouthWest,
        __SouthEast.name() : __SouthEast,
        __Center.name() : __Center
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.GroundCornersPointsType = GroundCornersPointsType
Namespace.addCategoryObject('typeBinding', 'GroundCornersPointsType', GroundCornersPointsType)


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_103 (pyxb.binding.basis.complexTypeDefinition):
    """Complex type [anonymous] with content type SIMPLE"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1590, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_103_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1593, 18)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1593, 18)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_103 = CTD_ANON_103


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_104 (pyxb.binding.basis.complexTypeDefinition):
    """Complex type [anonymous] with content type SIMPLE"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1599, 12)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_104_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1602, 18)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1602, 18)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_104 = CTD_ANON_104


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_105 (pyxb.binding.basis.complexTypeDefinition):
    """Complex type [anonymous] with content type SIMPLE"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1649, 8)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.double
    
    # Attribute unit uses Python identifier unit
    __unit = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'unit'), 'unit', '__aresysTypes_CTD_ANON_105_unit', _module_typeBindings.units)
    __unit._DeclarationLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1652, 14)
    __unit._UseLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1652, 14)
    
    unit = property(__unit.value, __unit.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __unit.name() : __unit
    })
_module_typeBindings.CTD_ANON_105 = CTD_ANON_105


# Complex type {aresysTypes}BurstInfoType with content type ELEMENT_ONLY
class BurstInfoType (TreeElementBaseType):
    """Bursts information"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'BurstInfoType')
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1659, 2)
    _ElementMap = TreeElementBaseType._ElementMap.copy()
    _AttributeMap = TreeElementBaseType._AttributeMap.copy()
    # Base type is TreeElementBaseType
    
    # Element NumberOfBursts uses Python identifier NumberOfBursts
    __NumberOfBursts = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'NumberOfBursts'), 'NumberOfBursts', '__aresysTypes_BurstInfoType_NumberOfBursts', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1666, 10), )

    
    NumberOfBursts = property(__NumberOfBursts.value, __NumberOfBursts.set, None, 'Number of bursts in the swath')

    
    # Element LinesPerBurst uses Python identifier LinesPerBurst
    __LinesPerBurst = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'LinesPerBurst'), 'LinesPerBurst', '__aresysTypes_BurstInfoType_LinesPerBurst', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1672, 12), )

    
    LinesPerBurst = property(__LinesPerBurst.value, __LinesPerBurst.set, None, 'Number of lines in each burst')

    
    # Element LinesPerBurstChangeList uses Python identifier LinesPerBurstChangeList
    __LinesPerBurstChangeList = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'LinesPerBurstChangeList'), 'LinesPerBurstChangeList', '__aresysTypes_BurstInfoType_LinesPerBurstChangeList', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1677, 12), )

    
    LinesPerBurstChangeList = property(__LinesPerBurstChangeList.value, __LinesPerBurstChangeList.set, None, None)

    
    # Element BurstRepetitionFrequency uses Python identifier BurstRepetitionFrequency
    __BurstRepetitionFrequency = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'BurstRepetitionFrequency'), 'BurstRepetitionFrequency', '__aresysTypes_BurstInfoType_BurstRepetitionFrequency', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1693, 10), )

    
    BurstRepetitionFrequency = property(__BurstRepetitionFrequency.value, __BurstRepetitionFrequency.set, None, 'Burst repetition frequency [Hz]')

    
    # Element Burst uses Python identifier Burst
    __Burst = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'Burst'), 'Burst', '__aresysTypes_BurstInfoType_Burst', True, pyxb.utils.utility.Location('aresysTypes.xsd', 1698, 10), )

    
    Burst = property(__Burst.value, __Burst.set, None, 'Time coordinates of each burst')

    
    # Attribute Number inherited from {aresysTypes}TreeElementBaseType
    
    # Attribute Total inherited from {aresysTypes}TreeElementBaseType
    _ElementMap.update({
        __NumberOfBursts.name() : __NumberOfBursts,
        __LinesPerBurst.name() : __LinesPerBurst,
        __LinesPerBurstChangeList.name() : __LinesPerBurstChangeList,
        __BurstRepetitionFrequency.name() : __BurstRepetitionFrequency,
        __Burst.name() : __Burst
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.BurstInfoType = BurstInfoType
Namespace.addCategoryObject('typeBinding', 'BurstInfoType', BurstInfoType)


# Complex type {aresysTypes}AntennaInfoType with content type ELEMENT_ONLY
class AntennaInfoType (TreeElementBaseType):
    """Antenna pattern information"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'AntennaInfoType')
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 1707, 2)
    _ElementMap = TreeElementBaseType._ElementMap.copy()
    _AttributeMap = TreeElementBaseType._AttributeMap.copy()
    # Base type is TreeElementBaseType
    
    # Element SensorName uses Python identifier SensorName
    __SensorName = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'SensorName'), 'SensorName', '__aresysTypes_AntennaInfoType_SensorName', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1714, 10), )

    
    SensorName = property(__SensorName.value, __SensorName.set, None, 'Sensor name: ASAR, PALSAR, ...')

    
    # Element AcquisitionMode uses Python identifier AcquisitionMode
    __AcquisitionMode = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'AcquisitionMode'), 'AcquisitionMode', '__aresysTypes_AntennaInfoType_AcquisitionMode', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1719, 10), )

    
    AcquisitionMode = property(__AcquisitionMode.value, __AcquisitionMode.set, None, 'Acquisition mode: STRIPMAP, TOPSAR, ...')

    
    # Element BeamName uses Python identifier BeamName
    __BeamName = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'BeamName'), 'BeamName', '__aresysTypes_AntennaInfoType_BeamName', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1724, 10), )

    
    BeamName = property(__BeamName.value, __BeamName.set, None, 'Acquisition beam name')

    
    # Element Polarization uses Python identifier Polarization
    __Polarization = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'Polarization'), 'Polarization', '__aresysTypes_AntennaInfoType_Polarization', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1729, 10), )

    
    Polarization = property(__Polarization.value, __Polarization.set, None, 'Antenna polarization (H/H, H/V...)')

    
    # Element LinesPerPattern uses Python identifier LinesPerPattern
    __LinesPerPattern = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'LinesPerPattern'), 'LinesPerPattern', '__aresysTypes_AntennaInfoType_LinesPerPattern', False, pyxb.utils.utility.Location('aresysTypes.xsd', 1734, 10), )

    
    LinesPerPattern = property(__LinesPerPattern.value, __LinesPerPattern.set, None, 'Contains number of lines for each pattern')

    
    # Attribute Number inherited from {aresysTypes}TreeElementBaseType
    
    # Attribute Total inherited from {aresysTypes}TreeElementBaseType
    _ElementMap.update({
        __SensorName.name() : __SensorName,
        __AcquisitionMode.name() : __AcquisitionMode,
        __BeamName.name() : __BeamName,
        __Polarization.name() : __Polarization,
        __LinesPerPattern.name() : __LinesPerPattern
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.AntennaInfoType = AntennaInfoType
Namespace.addCategoryObject('typeBinding', 'AntennaInfoType', AntennaInfoType)


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_106 (doubleWithUnit):
    """Complex type [anonymous] with content type SIMPLE"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 185, 12)
    _ElementMap = doubleWithUnit._ElementMap.copy()
    _AttributeMap = doubleWithUnit._AttributeMap.copy()
    # Base type is doubleWithUnit
    
    # Attribute unit inherited from {aresysTypes}doubleWithUnit
    _ElementMap.update({
        
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.CTD_ANON_106 = CTD_ANON_106


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_107 (doubleWithUnit):
    """Complex type [anonymous] with content type SIMPLE"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 192, 12)
    _ElementMap = doubleWithUnit._ElementMap.copy()
    _AttributeMap = doubleWithUnit._AttributeMap.copy()
    # Base type is doubleWithUnit
    
    # Attribute unit inherited from {aresysTypes}doubleWithUnit
    _ElementMap.update({
        
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.CTD_ANON_107 = CTD_ANON_107


# Complex type [anonymous] with content type SIMPLE
class CTD_ANON_108 (doubleWithUnit):
    """Complex type [anonymous] with content type SIMPLE"""
    _TypeDefinition = pyxb.binding.datatypes.double
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresysTypes.xsd', 199, 12)
    _ElementMap = doubleWithUnit._ElementMap.copy()
    _AttributeMap = doubleWithUnit._AttributeMap.copy()
    # Base type is doubleWithUnit
    
    # Attribute unit inherited from {aresysTypes}doubleWithUnit
    _ElementMap.update({
        
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.CTD_ANON_108 = CTD_ANON_108




CTD_ANON_3._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'val'), CTD_ANON_46, scope=CTD_ANON_3, location=pyxb.utils.utility.Location('aresysTypes.xsd', 472, 16)))

def _BuildAutomaton ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton
    del _BuildAutomaton
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=0, max=None, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 472, 16))
    counters.add(cc_0)
    states = []
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_0, False))
    symbol = pyxb.binding.content.ElementUse(CTD_ANON_3._UseForTag(pyxb.namespace.ExpandedName(None, 'val')), pyxb.utils.utility.Location('aresysTypes.xsd', 472, 16))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    transitions.append(fac.Transition(st_0, [
        fac.UpdateInstruction(cc_0, True) ]))
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, True, containing_state=None)
CTD_ANON_3._Automaton = _BuildAutomaton()




CTD_ANON_4._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'val'), CTD_ANON_47, scope=CTD_ANON_4, location=pyxb.utils.utility.Location('aresysTypes.xsd', 491, 16)))

def _BuildAutomaton_ ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_
    del _BuildAutomaton_
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=0, max=None, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 491, 16))
    counters.add(cc_0)
    states = []
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_0, False))
    symbol = pyxb.binding.content.ElementUse(CTD_ANON_4._UseForTag(pyxb.namespace.ExpandedName(None, 'val')), pyxb.utils.utility.Location('aresysTypes.xsd', 491, 16))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    transitions.append(fac.Transition(st_0, [
        fac.UpdateInstruction(cc_0, True) ]))
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, True, containing_state=None)
CTD_ANON_4._Automaton = _BuildAutomaton_()




CTD_ANON_5._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'val'), CTD_ANON_50, scope=CTD_ANON_5, location=pyxb.utils.utility.Location('aresysTypes.xsd', 544, 16)))

def _BuildAutomaton_2 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_2
    del _BuildAutomaton_2
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=3, max=3, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 544, 16))
    counters.add(cc_0)
    states = []
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_0, False))
    symbol = pyxb.binding.content.ElementUse(CTD_ANON_5._UseForTag(pyxb.namespace.ExpandedName(None, 'val')), pyxb.utils.utility.Location('aresysTypes.xsd', 544, 16))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    transitions.append(fac.Transition(st_0, [
        fac.UpdateInstruction(cc_0, True) ]))
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=None)
CTD_ANON_5._Automaton = _BuildAutomaton_2()




CTD_ANON_6._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'val'), CTD_ANON_53, scope=CTD_ANON_6, location=pyxb.utils.utility.Location('aresysTypes.xsd', 603, 16)))

def _BuildAutomaton_3 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_3
    del _BuildAutomaton_3
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=0, max=None, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 603, 16))
    counters.add(cc_0)
    states = []
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_0, False))
    symbol = pyxb.binding.content.ElementUse(CTD_ANON_6._UseForTag(pyxb.namespace.ExpandedName(None, 'val')), pyxb.utils.utility.Location('aresysTypes.xsd', 603, 16))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    transitions.append(fac.Transition(st_0, [
        fac.UpdateInstruction(cc_0, True) ]))
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, True, containing_state=None)
CTD_ANON_6._Automaton = _BuildAutomaton_3()




CTD_ANON_7._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'val'), CTD_ANON_54, scope=CTD_ANON_7, location=pyxb.utils.utility.Location('aresysTypes.xsd', 622, 16)))

def _BuildAutomaton_4 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_4
    del _BuildAutomaton_4
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=0, max=None, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 622, 16))
    counters.add(cc_0)
    states = []
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_0, False))
    symbol = pyxb.binding.content.ElementUse(CTD_ANON_7._UseForTag(pyxb.namespace.ExpandedName(None, 'val')), pyxb.utils.utility.Location('aresysTypes.xsd', 622, 16))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    transitions.append(fac.Transition(st_0, [
        fac.UpdateInstruction(cc_0, True) ]))
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, True, containing_state=None)
CTD_ANON_7._Automaton = _BuildAutomaton_4()




CTD_ANON_8._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'val'), CTD_ANON_55, scope=CTD_ANON_8, location=pyxb.utils.utility.Location('aresysTypes.xsd', 641, 16)))

def _BuildAutomaton_5 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_5
    del _BuildAutomaton_5
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=0, max=None, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 641, 16))
    counters.add(cc_0)
    states = []
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_0, False))
    symbol = pyxb.binding.content.ElementUse(CTD_ANON_8._UseForTag(pyxb.namespace.ExpandedName(None, 'val')), pyxb.utils.utility.Location('aresysTypes.xsd', 641, 16))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    transitions.append(fac.Transition(st_0, [
        fac.UpdateInstruction(cc_0, True) ]))
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, True, containing_state=None)
CTD_ANON_8._Automaton = _BuildAutomaton_5()




CTD_ANON_9._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'val'), CTD_ANON_61, scope=CTD_ANON_9, location=pyxb.utils.utility.Location('aresysTypes.xsd', 756, 16)))

def _BuildAutomaton_6 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_6
    del _BuildAutomaton_6
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=3, max=3, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 756, 16))
    counters.add(cc_0)
    states = []
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_0, False))
    symbol = pyxb.binding.content.ElementUse(CTD_ANON_9._UseForTag(pyxb.namespace.ExpandedName(None, 'val')), pyxb.utils.utility.Location('aresysTypes.xsd', 756, 16))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    transitions.append(fac.Transition(st_0, [
        fac.UpdateInstruction(cc_0, True) ]))
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=None)
CTD_ANON_9._Automaton = _BuildAutomaton_6()




CTD_ANON_11._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'val'), CTD_ANON_62, scope=CTD_ANON_11, location=pyxb.utils.utility.Location('aresysTypes.xsd', 816, 16)))

def _BuildAutomaton_7 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_7
    del _BuildAutomaton_7
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=0, max=None, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 816, 16))
    counters.add(cc_0)
    states = []
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_0, False))
    symbol = pyxb.binding.content.ElementUse(CTD_ANON_11._UseForTag(pyxb.namespace.ExpandedName(None, 'val')), pyxb.utils.utility.Location('aresysTypes.xsd', 816, 16))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    transitions.append(fac.Transition(st_0, [
        fac.UpdateInstruction(cc_0, True) ]))
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, True, containing_state=None)
CTD_ANON_11._Automaton = _BuildAutomaton_7()




CTD_ANON_13._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'val'), CTD_ANON_63, scope=CTD_ANON_13, location=pyxb.utils.utility.Location('aresysTypes.xsd', 845, 16)))

def _BuildAutomaton_8 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_8
    del _BuildAutomaton_8
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=0, max=None, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 845, 16))
    counters.add(cc_0)
    states = []
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_0, False))
    symbol = pyxb.binding.content.ElementUse(CTD_ANON_13._UseForTag(pyxb.namespace.ExpandedName(None, 'val')), pyxb.utils.utility.Location('aresysTypes.xsd', 845, 16))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    transitions.append(fac.Transition(st_0, [
        fac.UpdateInstruction(cc_0, True) ]))
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, True, containing_state=None)
CTD_ANON_13._Automaton = _BuildAutomaton_8()




CTD_ANON_15._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'val'), CTD_ANON_64, scope=CTD_ANON_15, location=pyxb.utils.utility.Location('aresysTypes.xsd', 874, 22)))

def _BuildAutomaton_9 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_9
    del _BuildAutomaton_9
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=0, max=None, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 874, 22))
    counters.add(cc_0)
    states = []
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_0, False))
    symbol = pyxb.binding.content.ElementUse(CTD_ANON_15._UseForTag(pyxb.namespace.ExpandedName(None, 'val')), pyxb.utils.utility.Location('aresysTypes.xsd', 874, 22))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    transitions.append(fac.Transition(st_0, [
        fac.UpdateInstruction(cc_0, True) ]))
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, True, containing_state=None)
CTD_ANON_15._Automaton = _BuildAutomaton_9()




CTD_ANON_16._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'val'), CTD_ANON_65, scope=CTD_ANON_16, location=pyxb.utils.utility.Location('aresysTypes.xsd', 893, 22)))

def _BuildAutomaton_10 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_10
    del _BuildAutomaton_10
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=0, max=None, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 893, 22))
    counters.add(cc_0)
    states = []
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_0, False))
    symbol = pyxb.binding.content.ElementUse(CTD_ANON_16._UseForTag(pyxb.namespace.ExpandedName(None, 'val')), pyxb.utils.utility.Location('aresysTypes.xsd', 893, 22))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    transitions.append(fac.Transition(st_0, [
        fac.UpdateInstruction(cc_0, True) ]))
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, True, containing_state=None)
CTD_ANON_16._Automaton = _BuildAutomaton_10()




CTD_ANON_18._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'val'), CTD_ANON_66, scope=CTD_ANON_18, location=pyxb.utils.utility.Location('aresysTypes.xsd', 922, 16)))

def _BuildAutomaton_11 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_11
    del _BuildAutomaton_11
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=0, max=None, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 922, 16))
    counters.add(cc_0)
    states = []
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_0, False))
    symbol = pyxb.binding.content.ElementUse(CTD_ANON_18._UseForTag(pyxb.namespace.ExpandedName(None, 'val')), pyxb.utils.utility.Location('aresysTypes.xsd', 922, 16))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    transitions.append(fac.Transition(st_0, [
        fac.UpdateInstruction(cc_0, True) ]))
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, True, containing_state=None)
CTD_ANON_18._Automaton = _BuildAutomaton_11()




CTD_ANON_19._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'val'), CTD_ANON_67, scope=CTD_ANON_19, location=pyxb.utils.utility.Location('aresysTypes.xsd', 941, 16)))

def _BuildAutomaton_12 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_12
    del _BuildAutomaton_12
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=0, max=None, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 941, 16))
    counters.add(cc_0)
    states = []
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_0, False))
    symbol = pyxb.binding.content.ElementUse(CTD_ANON_19._UseForTag(pyxb.namespace.ExpandedName(None, 'val')), pyxb.utils.utility.Location('aresysTypes.xsd', 941, 16))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    transitions.append(fac.Transition(st_0, [
        fac.UpdateInstruction(cc_0, True) ]))
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, True, containing_state=None)
CTD_ANON_19._Automaton = _BuildAutomaton_12()




CTD_ANON_20._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'val'), CTD_ANON_68, scope=CTD_ANON_20, location=pyxb.utils.utility.Location('aresysTypes.xsd', 965, 16)))

def _BuildAutomaton_13 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_13
    del _BuildAutomaton_13
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=0, max=None, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 965, 16))
    counters.add(cc_0)
    states = []
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_0, False))
    symbol = pyxb.binding.content.ElementUse(CTD_ANON_20._UseForTag(pyxb.namespace.ExpandedName(None, 'val')), pyxb.utils.utility.Location('aresysTypes.xsd', 965, 16))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    transitions.append(fac.Transition(st_0, [
        fac.UpdateInstruction(cc_0, True) ]))
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, True, containing_state=None)
CTD_ANON_20._Automaton = _BuildAutomaton_13()




CTD_ANON_21._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'val'), CTD_ANON_69, scope=CTD_ANON_21, location=pyxb.utils.utility.Location('aresysTypes.xsd', 989, 16)))

def _BuildAutomaton_14 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_14
    del _BuildAutomaton_14
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=0, max=None, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 989, 16))
    counters.add(cc_0)
    states = []
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_0, False))
    symbol = pyxb.binding.content.ElementUse(CTD_ANON_21._UseForTag(pyxb.namespace.ExpandedName(None, 'val')), pyxb.utils.utility.Location('aresysTypes.xsd', 989, 16))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    transitions.append(fac.Transition(st_0, [
        fac.UpdateInstruction(cc_0, True) ]))
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, True, containing_state=None)
CTD_ANON_21._Automaton = _BuildAutomaton_14()




CTD_ANON_22._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'val'), CTD_ANON_70, scope=CTD_ANON_22, location=pyxb.utils.utility.Location('aresysTypes.xsd', 1014, 18)))

def _BuildAutomaton_15 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_15
    del _BuildAutomaton_15
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=0, max=None, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 1014, 18))
    counters.add(cc_0)
    states = []
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_0, False))
    symbol = pyxb.binding.content.ElementUse(CTD_ANON_22._UseForTag(pyxb.namespace.ExpandedName(None, 'val')), pyxb.utils.utility.Location('aresysTypes.xsd', 1014, 18))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    transitions.append(fac.Transition(st_0, [
        fac.UpdateInstruction(cc_0, True) ]))
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, True, containing_state=None)
CTD_ANON_22._Automaton = _BuildAutomaton_15()




CTD_ANON_23._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'val'), CTD_ANON_71, scope=CTD_ANON_23, location=pyxb.utils.utility.Location('aresysTypes.xsd', 1033, 18)))

def _BuildAutomaton_16 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_16
    del _BuildAutomaton_16
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=0, max=None, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 1033, 18))
    counters.add(cc_0)
    states = []
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_0, False))
    symbol = pyxb.binding.content.ElementUse(CTD_ANON_23._UseForTag(pyxb.namespace.ExpandedName(None, 'val')), pyxb.utils.utility.Location('aresysTypes.xsd', 1033, 18))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    transitions.append(fac.Transition(st_0, [
        fac.UpdateInstruction(cc_0, True) ]))
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, True, containing_state=None)
CTD_ANON_23._Automaton = _BuildAutomaton_16()




CTD_ANON_24._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'DataBlockStatistic'), DataBlockStatisticsType, scope=CTD_ANON_24, location=pyxb.utils.utility.Location('aresysTypes.xsd', 1329, 16)))

def _BuildAutomaton_17 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_17
    del _BuildAutomaton_17
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(CTD_ANON_24._UseForTag(pyxb.namespace.ExpandedName(None, 'DataBlockStatistic')), pyxb.utils.utility.Location('aresysTypes.xsd', 1329, 16))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    transitions.append(fac.Transition(st_0, [
         ]))
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=None)
CTD_ANON_24._Automaton = _BuildAutomaton_17()




CTD_ANON_25._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'val'), CTD_ANON_96, scope=CTD_ANON_25, location=pyxb.utils.utility.Location('aresysTypes.xsd', 1471, 16)))

def _BuildAutomaton_18 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_18
    del _BuildAutomaton_18
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=7, max=None, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 1471, 16))
    counters.add(cc_0)
    states = []
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_0, False))
    symbol = pyxb.binding.content.ElementUse(CTD_ANON_25._UseForTag(pyxb.namespace.ExpandedName(None, 'val')), pyxb.utils.utility.Location('aresysTypes.xsd', 1471, 16))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    transitions.append(fac.Transition(st_0, [
        fac.UpdateInstruction(cc_0, True) ]))
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=None)
CTD_ANON_25._Automaton = _BuildAutomaton_18()




CTD_ANON_26._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'val'), CTD_ANON_99, scope=CTD_ANON_26, location=pyxb.utils.utility.Location('aresysTypes.xsd', 1525, 16)))

def _BuildAutomaton_19 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_19
    del _BuildAutomaton_19
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=4, max=4, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 1525, 16))
    counters.add(cc_0)
    states = []
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_0, False))
    symbol = pyxb.binding.content.ElementUse(CTD_ANON_26._UseForTag(pyxb.namespace.ExpandedName(None, 'val')), pyxb.utils.utility.Location('aresysTypes.xsd', 1525, 16))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    transitions.append(fac.Transition(st_0, [
        fac.UpdateInstruction(cc_0, True) ]))
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=None)
CTD_ANON_26._Automaton = _BuildAutomaton_19()




CTD_ANON_27._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'val'), CTD_ANON_100, scope=CTD_ANON_27, location=pyxb.utils.utility.Location('aresysTypes.xsd', 1544, 16)))

def _BuildAutomaton_20 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_20
    del _BuildAutomaton_20
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=4, max=4, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 1544, 16))
    counters.add(cc_0)
    states = []
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_0, False))
    symbol = pyxb.binding.content.ElementUse(CTD_ANON_27._UseForTag(pyxb.namespace.ExpandedName(None, 'val')), pyxb.utils.utility.Location('aresysTypes.xsd', 1544, 16))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    transitions.append(fac.Transition(st_0, [
        fac.UpdateInstruction(cc_0, True) ]))
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=None)
CTD_ANON_27._Automaton = _BuildAutomaton_20()




CTD_ANON_28._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'Point'), PointType, scope=CTD_ANON_28, location=pyxb.utils.utility.Location('aresysTypes.xsd', 1610, 16)))

def _BuildAutomaton_21 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_21
    del _BuildAutomaton_21
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(CTD_ANON_28._UseForTag(pyxb.namespace.ExpandedName(None, 'Point')), pyxb.utils.utility.Location('aresysTypes.xsd', 1610, 16))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=None)
CTD_ANON_28._Automaton = _BuildAutomaton_21()




CTD_ANON_29._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'Point'), PointType, scope=CTD_ANON_29, location=pyxb.utils.utility.Location('aresysTypes.xsd', 1617, 16)))

def _BuildAutomaton_22 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_22
    del _BuildAutomaton_22
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(CTD_ANON_29._UseForTag(pyxb.namespace.ExpandedName(None, 'Point')), pyxb.utils.utility.Location('aresysTypes.xsd', 1617, 16))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=None)
CTD_ANON_29._Automaton = _BuildAutomaton_22()




CTD_ANON_30._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'Point'), PointType, scope=CTD_ANON_30, location=pyxb.utils.utility.Location('aresysTypes.xsd', 1624, 16)))

def _BuildAutomaton_23 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_23
    del _BuildAutomaton_23
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(CTD_ANON_30._UseForTag(pyxb.namespace.ExpandedName(None, 'Point')), pyxb.utils.utility.Location('aresysTypes.xsd', 1624, 16))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=None)
CTD_ANON_30._Automaton = _BuildAutomaton_23()




CTD_ANON_31._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'Point'), PointType, scope=CTD_ANON_31, location=pyxb.utils.utility.Location('aresysTypes.xsd', 1631, 16)))

def _BuildAutomaton_24 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_24
    del _BuildAutomaton_24
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(CTD_ANON_31._UseForTag(pyxb.namespace.ExpandedName(None, 'Point')), pyxb.utils.utility.Location('aresysTypes.xsd', 1631, 16))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=None)
CTD_ANON_31._Automaton = _BuildAutomaton_24()




CTD_ANON_32._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'Point'), PointType, scope=CTD_ANON_32, location=pyxb.utils.utility.Location('aresysTypes.xsd', 1638, 16)))

def _BuildAutomaton_25 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_25
    del _BuildAutomaton_25
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(CTD_ANON_32._UseForTag(pyxb.namespace.ExpandedName(None, 'Point')), pyxb.utils.utility.Location('aresysTypes.xsd', 1638, 16))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=None)
CTD_ANON_32._Automaton = _BuildAutomaton_25()




PointType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'val'), CTD_ANON_105, scope=PointType, location=pyxb.utils.utility.Location('aresysTypes.xsd', 1648, 6)))

def _BuildAutomaton_26 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_26
    del _BuildAutomaton_26
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=5, max=5, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 1647, 4))
    counters.add(cc_0)
    states = []
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_0, False))
    symbol = pyxb.binding.content.ElementUse(PointType._UseForTag(pyxb.namespace.ExpandedName(None, 'val')), pyxb.utils.utility.Location('aresysTypes.xsd', 1648, 6))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    transitions.append(fac.Transition(st_0, [
        fac.UpdateInstruction(cc_0, True) ]))
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=None)
PointType._Automaton = _BuildAutomaton_26()




CTD_ANON_33._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'Lines'), CTD_ANON_34, scope=CTD_ANON_33, location=pyxb.utils.utility.Location('aresysTypes.xsd', 1680, 24)))

def _BuildAutomaton_27 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_27
    del _BuildAutomaton_27
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(CTD_ANON_33._UseForTag(pyxb.namespace.ExpandedName(None, 'Lines')), pyxb.utils.utility.Location('aresysTypes.xsd', 1680, 24))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    transitions.append(fac.Transition(st_0, [
         ]))
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=None)
CTD_ANON_33._Automaton = _BuildAutomaton_27()




BurstType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'RangeStartTime'), doubleWithUnit, scope=BurstType, documentation='Range absolute start time [s]', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1745, 6)))

BurstType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'AzimuthStartTime'), stringWithUnit, scope=BurstType, documentation='Azimuth start time absolute value [Utc]', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1750, 6)))

BurstType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'BurstCenterAzimuthShift'), doubleWithUnit, scope=BurstType, location=pyxb.utils.utility.Location('aresysTypes.xsd', 1755, 6)))

def _BuildAutomaton_28 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_28
    del _BuildAutomaton_28
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 1755, 6))
    counters.add(cc_0)
    states = []
    final_update = None
    symbol = pyxb.binding.content.ElementUse(BurstType._UseForTag(pyxb.namespace.ExpandedName(None, 'RangeStartTime')), pyxb.utils.utility.Location('aresysTypes.xsd', 1745, 6))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(BurstType._UseForTag(pyxb.namespace.ExpandedName(None, 'AzimuthStartTime')), pyxb.utils.utility.Location('aresysTypes.xsd', 1750, 6))
    st_1 = fac.State(symbol, is_initial=False, final_update=final_update, is_unordered_catenation=False)
    states.append(st_1)
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_0, False))
    symbol = pyxb.binding.content.ElementUse(BurstType._UseForTag(pyxb.namespace.ExpandedName(None, 'BurstCenterAzimuthShift')), pyxb.utils.utility.Location('aresysTypes.xsd', 1755, 6))
    st_2 = fac.State(symbol, is_initial=False, final_update=final_update, is_unordered_catenation=False)
    states.append(st_2)
    transitions = []
    transitions.append(fac.Transition(st_1, [
         ]))
    st_0._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_2, [
         ]))
    st_1._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_2, [
        fac.UpdateInstruction(cc_0, True) ]))
    st_2._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=None)
BurstType._Automaton = _BuildAutomaton_28()




DCOMPLEX._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'RealValue'), pyxb.binding.datatypes.double, scope=DCOMPLEX, location=pyxb.utils.utility.Location('aresysTypes.xsd', 1761, 6)))

DCOMPLEX._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'ImaginaryValue'), pyxb.binding.datatypes.double, scope=DCOMPLEX, location=pyxb.utils.utility.Location('aresysTypes.xsd', 1762, 6)))

def _BuildAutomaton_29 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_29
    del _BuildAutomaton_29
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = None
    symbol = pyxb.binding.content.ElementUse(DCOMPLEX._UseForTag(pyxb.namespace.ExpandedName(None, 'RealValue')), pyxb.utils.utility.Location('aresysTypes.xsd', 1761, 6))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(DCOMPLEX._UseForTag(pyxb.namespace.ExpandedName(None, 'ImaginaryValue')), pyxb.utils.utility.Location('aresysTypes.xsd', 1762, 6))
    st_1 = fac.State(symbol, is_initial=False, final_update=final_update, is_unordered_catenation=False)
    states.append(st_1)
    transitions = []
    transitions.append(fac.Transition(st_1, [
         ]))
    st_0._set_transitionSet(transitions)
    transitions = []
    st_1._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=None)
DCOMPLEX._Automaton = _BuildAutomaton_29()




ROIType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'SlowStartTime'), stringWithUnit, scope=ROIType, location=pyxb.utils.utility.Location('aresysTypes.xsd', 181, 10)))

ROIType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'SlowTimeLength'), doubleWithUnit, scope=ROIType, location=pyxb.utils.utility.Location('aresysTypes.xsd', 182, 10)))

ROIType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'SlowTimeDelta'), doubleWithUnit, scope=ROIType, location=pyxb.utils.utility.Location('aresysTypes.xsd', 183, 10)))

ROIType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'FastStartTime'), CTD_ANON_106, scope=ROIType, location=pyxb.utils.utility.Location('aresysTypes.xsd', 184, 10)))

ROIType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'FastTimeLength'), CTD_ANON_107, scope=ROIType, location=pyxb.utils.utility.Location('aresysTypes.xsd', 191, 10)))

ROIType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'FastTimeDelta'), CTD_ANON_108, scope=ROIType, location=pyxb.utils.utility.Location('aresysTypes.xsd', 198, 10)))

ROIType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'Lines'), CTD_ANON, scope=ROIType, location=pyxb.utils.utility.Location('aresysTypes.xsd', 205, 10)))

ROIType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'Samples'), CTD_ANON_, scope=ROIType, location=pyxb.utils.utility.Location('aresysTypes.xsd', 212, 10)))

def _BuildAutomaton_31 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_31
    del _BuildAutomaton_31
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(ROIType._UseForTag(pyxb.namespace.ExpandedName(None, 'SlowStartTime')), pyxb.utils.utility.Location('aresysTypes.xsd', 181, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_32 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_32
    del _BuildAutomaton_32
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(ROIType._UseForTag(pyxb.namespace.ExpandedName(None, 'SlowTimeLength')), pyxb.utils.utility.Location('aresysTypes.xsd', 182, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_33 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_33
    del _BuildAutomaton_33
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(ROIType._UseForTag(pyxb.namespace.ExpandedName(None, 'SlowTimeDelta')), pyxb.utils.utility.Location('aresysTypes.xsd', 183, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_34 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_34
    del _BuildAutomaton_34
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(ROIType._UseForTag(pyxb.namespace.ExpandedName(None, 'FastStartTime')), pyxb.utils.utility.Location('aresysTypes.xsd', 184, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_35 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_35
    del _BuildAutomaton_35
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(ROIType._UseForTag(pyxb.namespace.ExpandedName(None, 'FastTimeLength')), pyxb.utils.utility.Location('aresysTypes.xsd', 191, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_36 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_36
    del _BuildAutomaton_36
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(ROIType._UseForTag(pyxb.namespace.ExpandedName(None, 'FastTimeDelta')), pyxb.utils.utility.Location('aresysTypes.xsd', 198, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_37 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_37
    del _BuildAutomaton_37
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(ROIType._UseForTag(pyxb.namespace.ExpandedName(None, 'Lines')), pyxb.utils.utility.Location('aresysTypes.xsd', 205, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_38 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_38
    del _BuildAutomaton_38
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(ROIType._UseForTag(pyxb.namespace.ExpandedName(None, 'Samples')), pyxb.utils.utility.Location('aresysTypes.xsd', 212, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_30 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_30
    del _BuildAutomaton_30
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    sub_automata = []
    sub_automata.append(_BuildAutomaton_31())
    sub_automata.append(_BuildAutomaton_32())
    sub_automata.append(_BuildAutomaton_33())
    sub_automata.append(_BuildAutomaton_34())
    sub_automata.append(_BuildAutomaton_35())
    sub_automata.append(_BuildAutomaton_36())
    sub_automata.append(_BuildAutomaton_37())
    sub_automata.append(_BuildAutomaton_38())
    final_update = set()
    symbol = pyxb.utils.utility.Location('aresysTypes.xsd', 180, 8)
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=True)
    st_0._set_subAutomata(*sub_automata)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=None)
ROIType._Automaton = _BuildAutomaton_30()




PulseType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'Direction'), STD_ANON, scope=PulseType, documentation='Pulse direction (UP, DOWN)', location=pyxb.utils.utility.Location('aresysTypes.xsd', 230, 10)))

PulseType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'PulseLength'), doubleWithUnit, scope=PulseType, documentation='Pulse length [s]', location=pyxb.utils.utility.Location('aresysTypes.xsd', 241, 10)))

PulseType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'Bandwidth'), doubleWithUnit, scope=PulseType, documentation='Pulse bandwidth [Hz]', location=pyxb.utils.utility.Location('aresysTypes.xsd', 246, 10)))

PulseType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'PulseEnergy'), doubleWithUnit, scope=PulseType, documentation='Pulse energy [J]', location=pyxb.utils.utility.Location('aresysTypes.xsd', 251, 10)))

PulseType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'PulseSamplingRate'), doubleWithUnit, scope=PulseType, documentation='Pulse sampling rate [Hz]', location=pyxb.utils.utility.Location('aresysTypes.xsd', 256, 10)))

PulseType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'PulseStartFrequency'), doubleWithUnit, scope=PulseType, documentation='Pulse start frequency [Hz]', location=pyxb.utils.utility.Location('aresysTypes.xsd', 261, 10)))

PulseType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'PulseStartPhase'), doubleWithUnit, scope=PulseType, documentation='Pulse start phase [rad]', location=pyxb.utils.utility.Location('aresysTypes.xsd', 266, 10)))

def _BuildAutomaton_40 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_40
    del _BuildAutomaton_40
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 230, 10))
    counters.add(cc_0)
    states = []
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_0, False))
    symbol = pyxb.binding.content.ElementUse(PulseType._UseForTag(pyxb.namespace.ExpandedName(None, 'Direction')), pyxb.utils.utility.Location('aresysTypes.xsd', 230, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    transitions.append(fac.Transition(st_0, [
        fac.UpdateInstruction(cc_0, True) ]))
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, True, containing_state=st_0)

def _BuildAutomaton_41 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_41
    del _BuildAutomaton_41
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(PulseType._UseForTag(pyxb.namespace.ExpandedName(None, 'PulseLength')), pyxb.utils.utility.Location('aresysTypes.xsd', 241, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_42 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_42
    del _BuildAutomaton_42
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(PulseType._UseForTag(pyxb.namespace.ExpandedName(None, 'Bandwidth')), pyxb.utils.utility.Location('aresysTypes.xsd', 246, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_43 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_43
    del _BuildAutomaton_43
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(PulseType._UseForTag(pyxb.namespace.ExpandedName(None, 'PulseEnergy')), pyxb.utils.utility.Location('aresysTypes.xsd', 251, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_44 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_44
    del _BuildAutomaton_44
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(PulseType._UseForTag(pyxb.namespace.ExpandedName(None, 'PulseSamplingRate')), pyxb.utils.utility.Location('aresysTypes.xsd', 256, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_45 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_45
    del _BuildAutomaton_45
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 261, 10))
    counters.add(cc_0)
    states = []
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_0, False))
    symbol = pyxb.binding.content.ElementUse(PulseType._UseForTag(pyxb.namespace.ExpandedName(None, 'PulseStartFrequency')), pyxb.utils.utility.Location('aresysTypes.xsd', 261, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    transitions.append(fac.Transition(st_0, [
        fac.UpdateInstruction(cc_0, True) ]))
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, True, containing_state=st_0)

def _BuildAutomaton_46 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_46
    del _BuildAutomaton_46
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 266, 10))
    counters.add(cc_0)
    states = []
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_0, False))
    symbol = pyxb.binding.content.ElementUse(PulseType._UseForTag(pyxb.namespace.ExpandedName(None, 'PulseStartPhase')), pyxb.utils.utility.Location('aresysTypes.xsd', 266, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    transitions.append(fac.Transition(st_0, [
        fac.UpdateInstruction(cc_0, True) ]))
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, True, containing_state=st_0)

def _BuildAutomaton_39 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_39
    del _BuildAutomaton_39
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 230, 10))
    counters.add(cc_0)
    cc_1 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 261, 10))
    counters.add(cc_1)
    cc_2 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 266, 10))
    counters.add(cc_2)
    states = []
    sub_automata = []
    sub_automata.append(_BuildAutomaton_40())
    sub_automata.append(_BuildAutomaton_41())
    sub_automata.append(_BuildAutomaton_42())
    sub_automata.append(_BuildAutomaton_43())
    sub_automata.append(_BuildAutomaton_44())
    sub_automata.append(_BuildAutomaton_45())
    sub_automata.append(_BuildAutomaton_46())
    final_update = set()
    symbol = pyxb.utils.utility.Location('aresysTypes.xsd', 229, 8)
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=True)
    st_0._set_subAutomata(*sub_automata)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=None)
PulseType._Automaton = _BuildAutomaton_39()




DataSetInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'SensorName'), pyxb.binding.datatypes.string, scope=DataSetInfoType, documentation='Name of the sensor used to acquire the image: ASAR, PALSAR, ...', location=pyxb.utils.utility.Location('aresysTypes.xsd', 282, 10)))

DataSetInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'Description'), CTD_ANON_35, scope=DataSetInfoType, documentation='Description of the image', location=pyxb.utils.utility.Location('aresysTypes.xsd', 287, 10)))

DataSetInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'SenseDate'), CTD_ANON_36, scope=DataSetInfoType, documentation='Image acquisition date', location=pyxb.utils.utility.Location('aresysTypes.xsd', 299, 10)))

DataSetInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'AcquisitionMode'), CTD_ANON_37, scope=DataSetInfoType, documentation='Image acquisition mode: STRIPMAP, TOPSAR, ...', location=pyxb.utils.utility.Location('aresysTypes.xsd', 311, 10)))

DataSetInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'ImageType'), CTD_ANON_38, scope=DataSetInfoType, documentation='Image type: RAW DATA, RANGE FOCUSED, AZIMUTH FOCUSED ', location=pyxb.utils.utility.Location('aresysTypes.xsd', 323, 10)))

DataSetInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'Projection'), CTD_ANON_39, scope=DataSetInfoType, documentation='Image projection: SLANT RANGE, GROUND RANGE', location=pyxb.utils.utility.Location('aresysTypes.xsd', 335, 10)))

DataSetInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'ProjectionParameters'), CTD_ANON_2, scope=DataSetInfoType, location=pyxb.utils.utility.Location('aresysTypes.xsd', 347, 10)))

DataSetInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'AcquisitionStation'), CTD_ANON_40, scope=DataSetInfoType, documentation='Image acquisition station', location=pyxb.utils.utility.Location('aresysTypes.xsd', 356, 10)))

DataSetInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'ProcessingCenter'), CTD_ANON_41, scope=DataSetInfoType, documentation='Image processing center', location=pyxb.utils.utility.Location('aresysTypes.xsd', 368, 10)))

DataSetInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'ProcessingDate'), CTD_ANON_42, scope=DataSetInfoType, documentation='Image processing date', location=pyxb.utils.utility.Location('aresysTypes.xsd', 380, 10)))

DataSetInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'ProcessingSoftware'), CTD_ANON_43, scope=DataSetInfoType, documentation='Image processing software', location=pyxb.utils.utility.Location('aresysTypes.xsd', 392, 10)))

DataSetInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'fc_hz'), CTD_ANON_44, scope=DataSetInfoType, documentation='Radar carrier frequency [Hz]', location=pyxb.utils.utility.Location('aresysTypes.xsd', 404, 10)))

DataSetInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'SideLooking'), LeftRightType, scope=DataSetInfoType, documentation='Radar side looking: LEFT, RIGHT', location=pyxb.utils.utility.Location('aresysTypes.xsd', 416, 10)))

DataSetInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'ExternalCalibrationFactor'), pyxb.binding.datatypes.double, scope=DataSetInfoType, documentation='External calibration factor', location=pyxb.utils.utility.Location('aresysTypes.xsd', 421, 10)))

DataSetInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'DataTakeID'), pyxb.binding.datatypes.unsignedInt, scope=DataSetInfoType, location=pyxb.utils.utility.Location('aresysTypes.xsd', 426, 10)))

DataSetInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'InstrumentConfID'), STD_ANON_, scope=DataSetInfoType, location=pyxb.utils.utility.Location('aresysTypes.xsd', 427, 10)))

def _BuildAutomaton_48 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_48
    del _BuildAutomaton_48
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(DataSetInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'SensorName')), pyxb.utils.utility.Location('aresysTypes.xsd', 282, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_49 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_49
    del _BuildAutomaton_49
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(DataSetInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'Description')), pyxb.utils.utility.Location('aresysTypes.xsd', 287, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_50 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_50
    del _BuildAutomaton_50
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(DataSetInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'SenseDate')), pyxb.utils.utility.Location('aresysTypes.xsd', 299, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_51 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_51
    del _BuildAutomaton_51
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(DataSetInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'AcquisitionMode')), pyxb.utils.utility.Location('aresysTypes.xsd', 311, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_52 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_52
    del _BuildAutomaton_52
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(DataSetInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'ImageType')), pyxb.utils.utility.Location('aresysTypes.xsd', 323, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_53 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_53
    del _BuildAutomaton_53
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(DataSetInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'Projection')), pyxb.utils.utility.Location('aresysTypes.xsd', 335, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_54 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_54
    del _BuildAutomaton_54
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 347, 10))
    counters.add(cc_0)
    states = []
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_0, False))
    symbol = pyxb.binding.content.ElementUse(DataSetInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'ProjectionParameters')), pyxb.utils.utility.Location('aresysTypes.xsd', 347, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    transitions.append(fac.Transition(st_0, [
        fac.UpdateInstruction(cc_0, True) ]))
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, True, containing_state=st_0)

def _BuildAutomaton_55 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_55
    del _BuildAutomaton_55
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(DataSetInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'AcquisitionStation')), pyxb.utils.utility.Location('aresysTypes.xsd', 356, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_56 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_56
    del _BuildAutomaton_56
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(DataSetInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'ProcessingCenter')), pyxb.utils.utility.Location('aresysTypes.xsd', 368, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_57 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_57
    del _BuildAutomaton_57
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(DataSetInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'ProcessingDate')), pyxb.utils.utility.Location('aresysTypes.xsd', 380, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_58 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_58
    del _BuildAutomaton_58
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(DataSetInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'ProcessingSoftware')), pyxb.utils.utility.Location('aresysTypes.xsd', 392, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_59 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_59
    del _BuildAutomaton_59
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(DataSetInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'fc_hz')), pyxb.utils.utility.Location('aresysTypes.xsd', 404, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_60 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_60
    del _BuildAutomaton_60
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(DataSetInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'SideLooking')), pyxb.utils.utility.Location('aresysTypes.xsd', 416, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_61 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_61
    del _BuildAutomaton_61
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 421, 10))
    counters.add(cc_0)
    states = []
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_0, False))
    symbol = pyxb.binding.content.ElementUse(DataSetInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'ExternalCalibrationFactor')), pyxb.utils.utility.Location('aresysTypes.xsd', 421, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    transitions.append(fac.Transition(st_0, [
        fac.UpdateInstruction(cc_0, True) ]))
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, True, containing_state=st_0)

def _BuildAutomaton_62 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_62
    del _BuildAutomaton_62
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 426, 10))
    counters.add(cc_0)
    states = []
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_0, False))
    symbol = pyxb.binding.content.ElementUse(DataSetInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'DataTakeID')), pyxb.utils.utility.Location('aresysTypes.xsd', 426, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    transitions.append(fac.Transition(st_0, [
        fac.UpdateInstruction(cc_0, True) ]))
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, True, containing_state=st_0)

def _BuildAutomaton_63 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_63
    del _BuildAutomaton_63
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 427, 10))
    counters.add(cc_0)
    states = []
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_0, False))
    symbol = pyxb.binding.content.ElementUse(DataSetInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'InstrumentConfID')), pyxb.utils.utility.Location('aresysTypes.xsd', 427, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    transitions.append(fac.Transition(st_0, [
        fac.UpdateInstruction(cc_0, True) ]))
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, True, containing_state=st_0)

def _BuildAutomaton_47 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_47
    del _BuildAutomaton_47
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 347, 10))
    counters.add(cc_0)
    cc_1 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 421, 10))
    counters.add(cc_1)
    cc_2 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 426, 10))
    counters.add(cc_2)
    cc_3 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 427, 10))
    counters.add(cc_3)
    states = []
    sub_automata = []
    sub_automata.append(_BuildAutomaton_48())
    sub_automata.append(_BuildAutomaton_49())
    sub_automata.append(_BuildAutomaton_50())
    sub_automata.append(_BuildAutomaton_51())
    sub_automata.append(_BuildAutomaton_52())
    sub_automata.append(_BuildAutomaton_53())
    sub_automata.append(_BuildAutomaton_54())
    sub_automata.append(_BuildAutomaton_55())
    sub_automata.append(_BuildAutomaton_56())
    sub_automata.append(_BuildAutomaton_57())
    sub_automata.append(_BuildAutomaton_58())
    sub_automata.append(_BuildAutomaton_59())
    sub_automata.append(_BuildAutomaton_60())
    sub_automata.append(_BuildAutomaton_61())
    sub_automata.append(_BuildAutomaton_62())
    sub_automata.append(_BuildAutomaton_63())
    final_update = set()
    symbol = pyxb.utils.utility.Location('aresysTypes.xsd', 281, 8)
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=True)
    st_0._set_subAutomata(*sub_automata)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=None)
DataSetInfoType._Automaton = _BuildAutomaton_47()




StateVectorDataType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'OrbitNumber'), pyxb.binding.datatypes.string, scope=StateVectorDataType, documentation='Number of the orbit', location=pyxb.utils.utility.Location('aresysTypes.xsd', 446, 10)))

StateVectorDataType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'Track'), pyxb.binding.datatypes.string, scope=StateVectorDataType, documentation='Number of the track', location=pyxb.utils.utility.Location('aresysTypes.xsd', 451, 10)))

StateVectorDataType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'OrbitDirection'), CTD_ANON_45, scope=StateVectorDataType, documentation='Direction of the orbit: ASCENDING, DESCENDING', location=pyxb.utils.utility.Location('aresysTypes.xsd', 456, 10)))

StateVectorDataType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'pSV_m'), CTD_ANON_3, scope=StateVectorDataType, documentation='Orbit state vectors position coordinates (xyz) [m]', location=pyxb.utils.utility.Location('aresysTypes.xsd', 466, 10)))

StateVectorDataType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'vSV_mOs'), CTD_ANON_4, scope=StateVectorDataType, documentation='Orbit state vectors velocity coordinates [m/s]', location=pyxb.utils.utility.Location('aresysTypes.xsd', 485, 10)))

StateVectorDataType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 't_ref_Utc'), pyxb.binding.datatypes.string, scope=StateVectorDataType, documentation='Azimuth absolute start time for the first state vector [Utc]', location=pyxb.utils.utility.Location('aresysTypes.xsd', 504, 10)))

StateVectorDataType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'dtSV_s'), CTD_ANON_48, scope=StateVectorDataType, documentation='Azimuth time interval between two consecutive state vectors [s]', location=pyxb.utils.utility.Location('aresysTypes.xsd', 509, 10)))

StateVectorDataType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'nSV_n'), CTD_ANON_49, scope=StateVectorDataType, documentation='Number of state vectors', location=pyxb.utils.utility.Location('aresysTypes.xsd', 521, 10)))

StateVectorDataType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'AscendingNodeTime'), pyxb.binding.datatypes.string, scope=StateVectorDataType, documentation='Azimuth absolute time of the ascending node', location=pyxb.utils.utility.Location('aresysTypes.xsd', 533, 10)))

StateVectorDataType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'AscendingNodeCoords'), CTD_ANON_5, scope=StateVectorDataType, documentation='Coordinates of the ascending node', location=pyxb.utils.utility.Location('aresysTypes.xsd', 538, 10)))

def _BuildAutomaton_65 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_65
    del _BuildAutomaton_65
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(StateVectorDataType._UseForTag(pyxb.namespace.ExpandedName(None, 'OrbitNumber')), pyxb.utils.utility.Location('aresysTypes.xsd', 446, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_66 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_66
    del _BuildAutomaton_66
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(StateVectorDataType._UseForTag(pyxb.namespace.ExpandedName(None, 'Track')), pyxb.utils.utility.Location('aresysTypes.xsd', 451, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_67 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_67
    del _BuildAutomaton_67
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(StateVectorDataType._UseForTag(pyxb.namespace.ExpandedName(None, 'OrbitDirection')), pyxb.utils.utility.Location('aresysTypes.xsd', 456, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_68 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_68
    del _BuildAutomaton_68
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(StateVectorDataType._UseForTag(pyxb.namespace.ExpandedName(None, 'pSV_m')), pyxb.utils.utility.Location('aresysTypes.xsd', 466, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_69 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_69
    del _BuildAutomaton_69
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(StateVectorDataType._UseForTag(pyxb.namespace.ExpandedName(None, 'vSV_mOs')), pyxb.utils.utility.Location('aresysTypes.xsd', 485, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_70 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_70
    del _BuildAutomaton_70
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(StateVectorDataType._UseForTag(pyxb.namespace.ExpandedName(None, 't_ref_Utc')), pyxb.utils.utility.Location('aresysTypes.xsd', 504, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_71 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_71
    del _BuildAutomaton_71
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(StateVectorDataType._UseForTag(pyxb.namespace.ExpandedName(None, 'dtSV_s')), pyxb.utils.utility.Location('aresysTypes.xsd', 509, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_72 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_72
    del _BuildAutomaton_72
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(StateVectorDataType._UseForTag(pyxb.namespace.ExpandedName(None, 'nSV_n')), pyxb.utils.utility.Location('aresysTypes.xsd', 521, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_73 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_73
    del _BuildAutomaton_73
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 533, 10))
    counters.add(cc_0)
    states = []
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_0, False))
    symbol = pyxb.binding.content.ElementUse(StateVectorDataType._UseForTag(pyxb.namespace.ExpandedName(None, 'AscendingNodeTime')), pyxb.utils.utility.Location('aresysTypes.xsd', 533, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    transitions.append(fac.Transition(st_0, [
        fac.UpdateInstruction(cc_0, True) ]))
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, True, containing_state=st_0)

def _BuildAutomaton_74 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_74
    del _BuildAutomaton_74
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 538, 10))
    counters.add(cc_0)
    states = []
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_0, False))
    symbol = pyxb.binding.content.ElementUse(StateVectorDataType._UseForTag(pyxb.namespace.ExpandedName(None, 'AscendingNodeCoords')), pyxb.utils.utility.Location('aresysTypes.xsd', 538, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    transitions.append(fac.Transition(st_0, [
        fac.UpdateInstruction(cc_0, True) ]))
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, True, containing_state=st_0)

def _BuildAutomaton_64 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_64
    del _BuildAutomaton_64
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 533, 10))
    counters.add(cc_0)
    cc_1 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 538, 10))
    counters.add(cc_1)
    states = []
    sub_automata = []
    sub_automata.append(_BuildAutomaton_65())
    sub_automata.append(_BuildAutomaton_66())
    sub_automata.append(_BuildAutomaton_67())
    sub_automata.append(_BuildAutomaton_68())
    sub_automata.append(_BuildAutomaton_69())
    sub_automata.append(_BuildAutomaton_70())
    sub_automata.append(_BuildAutomaton_71())
    sub_automata.append(_BuildAutomaton_72())
    sub_automata.append(_BuildAutomaton_73())
    sub_automata.append(_BuildAutomaton_74())
    final_update = set()
    symbol = pyxb.utils.utility.Location('aresysTypes.xsd', 445, 8)
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=True)
    st_0._set_subAutomata(*sub_automata)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=None)
StateVectorDataType._Automaton = _BuildAutomaton_64()




AttitudeInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 't_ref_Utc'), pyxb.binding.datatypes.string, scope=AttitudeInfoType, documentation='Azimuth absolute start time for the first attitude value [Utc]', location=pyxb.utils.utility.Location('aresysTypes.xsd', 568, 10)))

AttitudeInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'dtYPR_s'), CTD_ANON_51, scope=AttitudeInfoType, documentation='Azimuth time interval between two consecutive attitude values [s]', location=pyxb.utils.utility.Location('aresysTypes.xsd', 573, 10)))

AttitudeInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'nYPR_n'), CTD_ANON_52, scope=AttitudeInfoType, documentation='Number of attitude values', location=pyxb.utils.utility.Location('aresysTypes.xsd', 585, 10)))

AttitudeInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'yaw_deg'), CTD_ANON_6, scope=AttitudeInfoType, documentation='Yaw angle values [deg]', location=pyxb.utils.utility.Location('aresysTypes.xsd', 597, 10)))

AttitudeInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'pitch_deg'), CTD_ANON_7, scope=AttitudeInfoType, documentation='Pitch angle values [deg]', location=pyxb.utils.utility.Location('aresysTypes.xsd', 616, 10)))

AttitudeInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'roll_deg'), CTD_ANON_8, scope=AttitudeInfoType, documentation='Roll angle values [deg]', location=pyxb.utils.utility.Location('aresysTypes.xsd', 635, 10)))

AttitudeInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'referenceFrame'), ReferenceFrameType, scope=AttitudeInfoType, documentation='Reference frame: GEOCENTRIC, GEODETIC, ZERODOPPLER', location=pyxb.utils.utility.Location('aresysTypes.xsd', 654, 10)))

AttitudeInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'rotationOrder'), RotationOrderType, scope=AttitudeInfoType, documentation='Rotation order: YPR, YRP, PRY, PYR, RPY, RYP', location=pyxb.utils.utility.Location('aresysTypes.xsd', 659, 10)))

AttitudeInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'AttitudeType'), AttitudeType, scope=AttitudeInfoType, documentation='Attitude type: NOMINAL, REFINED', location=pyxb.utils.utility.Location('aresysTypes.xsd', 664, 10)))

def _BuildAutomaton_76 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_76
    del _BuildAutomaton_76
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(AttitudeInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 't_ref_Utc')), pyxb.utils.utility.Location('aresysTypes.xsd', 568, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_77 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_77
    del _BuildAutomaton_77
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(AttitudeInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'dtYPR_s')), pyxb.utils.utility.Location('aresysTypes.xsd', 573, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_78 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_78
    del _BuildAutomaton_78
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(AttitudeInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'nYPR_n')), pyxb.utils.utility.Location('aresysTypes.xsd', 585, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_79 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_79
    del _BuildAutomaton_79
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(AttitudeInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'yaw_deg')), pyxb.utils.utility.Location('aresysTypes.xsd', 597, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_80 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_80
    del _BuildAutomaton_80
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(AttitudeInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'pitch_deg')), pyxb.utils.utility.Location('aresysTypes.xsd', 616, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_81 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_81
    del _BuildAutomaton_81
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(AttitudeInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'roll_deg')), pyxb.utils.utility.Location('aresysTypes.xsd', 635, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_82 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_82
    del _BuildAutomaton_82
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(AttitudeInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'referenceFrame')), pyxb.utils.utility.Location('aresysTypes.xsd', 654, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_83 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_83
    del _BuildAutomaton_83
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(AttitudeInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'rotationOrder')), pyxb.utils.utility.Location('aresysTypes.xsd', 659, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_84 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_84
    del _BuildAutomaton_84
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(AttitudeInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'AttitudeType')), pyxb.utils.utility.Location('aresysTypes.xsd', 664, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_75 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_75
    del _BuildAutomaton_75
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    sub_automata = []
    sub_automata.append(_BuildAutomaton_76())
    sub_automata.append(_BuildAutomaton_77())
    sub_automata.append(_BuildAutomaton_78())
    sub_automata.append(_BuildAutomaton_79())
    sub_automata.append(_BuildAutomaton_80())
    sub_automata.append(_BuildAutomaton_81())
    sub_automata.append(_BuildAutomaton_82())
    sub_automata.append(_BuildAutomaton_83())
    sub_automata.append(_BuildAutomaton_84())
    final_update = set()
    symbol = pyxb.utils.utility.Location('aresysTypes.xsd', 567, 8)
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=True)
    st_0._set_subAutomata(*sub_automata)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=None)
AttitudeInfoType._Automaton = _BuildAutomaton_75()




SwathInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'Swath'), CTD_ANON_56, scope=SwathInfoType, documentation='Swath name', location=pyxb.utils.utility.Location('aresysTypes.xsd', 680, 10)))

SwathInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'SwathAcquisitionOrder'), CTD_ANON_57, scope=SwathInfoType, documentation='Swath acquisition order', location=pyxb.utils.utility.Location('aresysTypes.xsd', 692, 10)))

SwathInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'Polarization'), PolarizationType, scope=SwathInfoType, documentation='Polarization: H/H, H/V, V/H, V/V', location=pyxb.utils.utility.Location('aresysTypes.xsd', 704, 10)))

SwathInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'Rank'), CTD_ANON_58, scope=SwathInfoType, documentation='Rank', location=pyxb.utils.utility.Location('aresysTypes.xsd', 709, 10)))

SwathInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'RangeDelayBias'), CTD_ANON_59, scope=SwathInfoType, documentation='Range delay bias [s]', location=pyxb.utils.utility.Location('aresysTypes.xsd', 721, 10)))

SwathInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'AcquisitionStartTime'), CTD_ANON_60, scope=SwathInfoType, documentation='Acquisition start time [Utc]', location=pyxb.utils.utility.Location('aresysTypes.xsd', 733, 10)))

SwathInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'AzimuthSteeringRateReferenceTime'), doubleWithUnit, scope=SwathInfoType, documentation='Azimuth antenna steering rate polynomial reference time [s]', location=pyxb.utils.utility.Location('aresysTypes.xsd', 745, 10)))

SwathInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'AzimuthSteeringRatePol'), CTD_ANON_9, scope=SwathInfoType, documentation='Azimuth antenna steering rate polynomial coefficients: const [rad/s], az [rad/s^2], az^2 [rad/s^3]', location=pyxb.utils.utility.Location('aresysTypes.xsd', 750, 10)))

SwathInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'AcquisitionPRF'), pyxb.binding.datatypes.double, scope=SwathInfoType, documentation='Acquisition Pulse Repetition Frequency', location=pyxb.utils.utility.Location('aresysTypes.xsd', 769, 10)))

SwathInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'EchoesPerBurst'), pyxb.binding.datatypes.unsignedInt, scope=SwathInfoType, documentation='Number of echoes for each burst', location=pyxb.utils.utility.Location('aresysTypes.xsd', 774, 10)))

SwathInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'ChannelDelay'), pyxb.binding.datatypes.double, scope=SwathInfoType, documentation='Range channel delay time', location=pyxb.utils.utility.Location('aresysTypes.xsd', 779, 10)))

SwathInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'RxGain'), pyxb.binding.datatypes.float, scope=SwathInfoType, documentation='Value of the commandable Rx attenuation in the receiver channel', location=pyxb.utils.utility.Location('aresysTypes.xsd', 784, 10)))

def _BuildAutomaton_86 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_86
    del _BuildAutomaton_86
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(SwathInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'Swath')), pyxb.utils.utility.Location('aresysTypes.xsd', 680, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_87 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_87
    del _BuildAutomaton_87
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(SwathInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'SwathAcquisitionOrder')), pyxb.utils.utility.Location('aresysTypes.xsd', 692, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_88 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_88
    del _BuildAutomaton_88
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(SwathInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'Polarization')), pyxb.utils.utility.Location('aresysTypes.xsd', 704, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_89 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_89
    del _BuildAutomaton_89
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(SwathInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'Rank')), pyxb.utils.utility.Location('aresysTypes.xsd', 709, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_90 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_90
    del _BuildAutomaton_90
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(SwathInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'RangeDelayBias')), pyxb.utils.utility.Location('aresysTypes.xsd', 721, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_91 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_91
    del _BuildAutomaton_91
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(SwathInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'AcquisitionStartTime')), pyxb.utils.utility.Location('aresysTypes.xsd', 733, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_92 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_92
    del _BuildAutomaton_92
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(SwathInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'AzimuthSteeringRateReferenceTime')), pyxb.utils.utility.Location('aresysTypes.xsd', 745, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_93 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_93
    del _BuildAutomaton_93
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(SwathInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'AzimuthSteeringRatePol')), pyxb.utils.utility.Location('aresysTypes.xsd', 750, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_94 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_94
    del _BuildAutomaton_94
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(SwathInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'AcquisitionPRF')), pyxb.utils.utility.Location('aresysTypes.xsd', 769, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_95 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_95
    del _BuildAutomaton_95
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(SwathInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'EchoesPerBurst')), pyxb.utils.utility.Location('aresysTypes.xsd', 774, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_96 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_96
    del _BuildAutomaton_96
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 779, 10))
    counters.add(cc_0)
    states = []
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_0, False))
    symbol = pyxb.binding.content.ElementUse(SwathInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'ChannelDelay')), pyxb.utils.utility.Location('aresysTypes.xsd', 779, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    transitions.append(fac.Transition(st_0, [
        fac.UpdateInstruction(cc_0, True) ]))
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, True, containing_state=st_0)

def _BuildAutomaton_97 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_97
    del _BuildAutomaton_97
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 784, 10))
    counters.add(cc_0)
    states = []
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_0, False))
    symbol = pyxb.binding.content.ElementUse(SwathInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'RxGain')), pyxb.utils.utility.Location('aresysTypes.xsd', 784, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    transitions.append(fac.Transition(st_0, [
        fac.UpdateInstruction(cc_0, True) ]))
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, True, containing_state=st_0)

def _BuildAutomaton_85 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_85
    del _BuildAutomaton_85
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 779, 10))
    counters.add(cc_0)
    cc_1 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 784, 10))
    counters.add(cc_1)
    states = []
    sub_automata = []
    sub_automata.append(_BuildAutomaton_86())
    sub_automata.append(_BuildAutomaton_87())
    sub_automata.append(_BuildAutomaton_88())
    sub_automata.append(_BuildAutomaton_89())
    sub_automata.append(_BuildAutomaton_90())
    sub_automata.append(_BuildAutomaton_91())
    sub_automata.append(_BuildAutomaton_92())
    sub_automata.append(_BuildAutomaton_93())
    sub_automata.append(_BuildAutomaton_94())
    sub_automata.append(_BuildAutomaton_95())
    sub_automata.append(_BuildAutomaton_96())
    sub_automata.append(_BuildAutomaton_97())
    final_update = set()
    symbol = pyxb.utils.utility.Location('aresysTypes.xsd', 679, 8)
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=True)
    st_0._set_subAutomata(*sub_automata)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=None)
SwathInfoType._Automaton = _BuildAutomaton_85()




AcquisitionTimelineType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'MissingLines_number'), CTD_ANON_10, scope=AcquisitionTimelineType, documentation='Number of missing lines', location=pyxb.utils.utility.Location('aresysTypes.xsd', 800, 10)))

AcquisitionTimelineType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'MissingLines_azimuthtimes'), CTD_ANON_11, scope=AcquisitionTimelineType, documentation='Azimuth relative times for each missing line', location=pyxb.utils.utility.Location('aresysTypes.xsd', 810, 10)))

AcquisitionTimelineType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'DuplicatedLines_number'), CTD_ANON_12, scope=AcquisitionTimelineType, documentation='Number of duplicated lines', location=pyxb.utils.utility.Location('aresysTypes.xsd', 829, 10)))

AcquisitionTimelineType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'DuplicatedLines_azimuthtimes'), CTD_ANON_13, scope=AcquisitionTimelineType, documentation='Azimuth relative times for each duplicated line', location=pyxb.utils.utility.Location('aresysTypes.xsd', 839, 10)))

AcquisitionTimelineType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'PRF_changes_number'), CTD_ANON_14, scope=AcquisitionTimelineType, documentation='Number of PRF changes', location=pyxb.utils.utility.Location('aresysTypes.xsd', 858, 10)))

AcquisitionTimelineType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'PRF_changes_azimuthtimes'), CTD_ANON_15, scope=AcquisitionTimelineType, documentation='Azimuth relative times for each PRF change', location=pyxb.utils.utility.Location('aresysTypes.xsd', 868, 10)))

AcquisitionTimelineType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'PRF_changes_values'), CTD_ANON_16, scope=AcquisitionTimelineType, documentation='PRF changes values', location=pyxb.utils.utility.Location('aresysTypes.xsd', 887, 10)))

AcquisitionTimelineType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'Swst_changes_number'), CTD_ANON_17, scope=AcquisitionTimelineType, documentation='Number of SWST changes', location=pyxb.utils.utility.Location('aresysTypes.xsd', 906, 10)))

AcquisitionTimelineType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'Swst_changes_azimuthtimes'), CTD_ANON_18, scope=AcquisitionTimelineType, documentation='Azimuth relative times for each SWST change', location=pyxb.utils.utility.Location('aresysTypes.xsd', 916, 10)))

AcquisitionTimelineType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'Swst_changes_values'), CTD_ANON_19, scope=AcquisitionTimelineType, documentation='SWST changes values', location=pyxb.utils.utility.Location('aresysTypes.xsd', 935, 10)))

AcquisitionTimelineType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'noise_packets_number'), pyxb.binding.datatypes.unsignedInt, scope=AcquisitionTimelineType, documentation='Number of noise packets', location=pyxb.utils.utility.Location('aresysTypes.xsd', 954, 10)))

AcquisitionTimelineType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'noise_packets_azimuthtimes'), CTD_ANON_20, scope=AcquisitionTimelineType, documentation='Azimuth relative times for each noise packet', location=pyxb.utils.utility.Location('aresysTypes.xsd', 959, 10)))

AcquisitionTimelineType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'Internal_calibration_number'), pyxb.binding.datatypes.unsignedInt, scope=AcquisitionTimelineType, documentation='Number of internal calibration packets', location=pyxb.utils.utility.Location('aresysTypes.xsd', 978, 10)))

AcquisitionTimelineType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'Internal_calibration_azimuthtimes'), CTD_ANON_21, scope=AcquisitionTimelineType, documentation='Azimuth relative times for each internal calibration packet', location=pyxb.utils.utility.Location('aresysTypes.xsd', 983, 10)))

AcquisitionTimelineType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'Swl_changes_number'), pyxb.binding.datatypes.unsignedInt, scope=AcquisitionTimelineType, documentation='Number of SWL changes', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1003, 12)))

AcquisitionTimelineType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'Swl_changes_azimuthtimes'), CTD_ANON_22, scope=AcquisitionTimelineType, documentation='Relative azimuth times for each SWL change', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1008, 12)))

AcquisitionTimelineType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'Swl_changes_values'), CTD_ANON_23, scope=AcquisitionTimelineType, documentation='SWL changes values', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1027, 12)))

AcquisitionTimelineType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'ChirpPeriod'), pyxb.binding.datatypes.string, scope=AcquisitionTimelineType, documentation='Periodic list of chirp indexes related to multi-chirp image acquisition', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1047, 10)))

def _BuildAutomaton_98 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_98
    del _BuildAutomaton_98
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 829, 10))
    counters.add(cc_0)
    cc_1 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 839, 10))
    counters.add(cc_1)
    cc_2 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 858, 10))
    counters.add(cc_2)
    cc_3 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 868, 10))
    counters.add(cc_3)
    cc_4 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 887, 10))
    counters.add(cc_4)
    cc_5 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 1002, 10))
    counters.add(cc_5)
    cc_6 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 1047, 10))
    counters.add(cc_6)
    states = []
    final_update = None
    symbol = pyxb.binding.content.ElementUse(AcquisitionTimelineType._UseForTag(pyxb.namespace.ExpandedName(None, 'MissingLines_number')), pyxb.utils.utility.Location('aresysTypes.xsd', 800, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    final_update = None
    symbol = pyxb.binding.content.ElementUse(AcquisitionTimelineType._UseForTag(pyxb.namespace.ExpandedName(None, 'MissingLines_azimuthtimes')), pyxb.utils.utility.Location('aresysTypes.xsd', 810, 10))
    st_1 = fac.State(symbol, is_initial=False, final_update=final_update, is_unordered_catenation=False)
    states.append(st_1)
    final_update = None
    symbol = pyxb.binding.content.ElementUse(AcquisitionTimelineType._UseForTag(pyxb.namespace.ExpandedName(None, 'DuplicatedLines_number')), pyxb.utils.utility.Location('aresysTypes.xsd', 829, 10))
    st_2 = fac.State(symbol, is_initial=False, final_update=final_update, is_unordered_catenation=False)
    states.append(st_2)
    final_update = None
    symbol = pyxb.binding.content.ElementUse(AcquisitionTimelineType._UseForTag(pyxb.namespace.ExpandedName(None, 'DuplicatedLines_azimuthtimes')), pyxb.utils.utility.Location('aresysTypes.xsd', 839, 10))
    st_3 = fac.State(symbol, is_initial=False, final_update=final_update, is_unordered_catenation=False)
    states.append(st_3)
    final_update = None
    symbol = pyxb.binding.content.ElementUse(AcquisitionTimelineType._UseForTag(pyxb.namespace.ExpandedName(None, 'PRF_changes_number')), pyxb.utils.utility.Location('aresysTypes.xsd', 858, 10))
    st_4 = fac.State(symbol, is_initial=False, final_update=final_update, is_unordered_catenation=False)
    states.append(st_4)
    final_update = None
    symbol = pyxb.binding.content.ElementUse(AcquisitionTimelineType._UseForTag(pyxb.namespace.ExpandedName(None, 'PRF_changes_azimuthtimes')), pyxb.utils.utility.Location('aresysTypes.xsd', 868, 10))
    st_5 = fac.State(symbol, is_initial=False, final_update=final_update, is_unordered_catenation=False)
    states.append(st_5)
    final_update = None
    symbol = pyxb.binding.content.ElementUse(AcquisitionTimelineType._UseForTag(pyxb.namespace.ExpandedName(None, 'PRF_changes_values')), pyxb.utils.utility.Location('aresysTypes.xsd', 887, 10))
    st_6 = fac.State(symbol, is_initial=False, final_update=final_update, is_unordered_catenation=False)
    states.append(st_6)
    final_update = None
    symbol = pyxb.binding.content.ElementUse(AcquisitionTimelineType._UseForTag(pyxb.namespace.ExpandedName(None, 'Swst_changes_number')), pyxb.utils.utility.Location('aresysTypes.xsd', 906, 10))
    st_7 = fac.State(symbol, is_initial=False, final_update=final_update, is_unordered_catenation=False)
    states.append(st_7)
    final_update = None
    symbol = pyxb.binding.content.ElementUse(AcquisitionTimelineType._UseForTag(pyxb.namespace.ExpandedName(None, 'Swst_changes_azimuthtimes')), pyxb.utils.utility.Location('aresysTypes.xsd', 916, 10))
    st_8 = fac.State(symbol, is_initial=False, final_update=final_update, is_unordered_catenation=False)
    states.append(st_8)
    final_update = None
    symbol = pyxb.binding.content.ElementUse(AcquisitionTimelineType._UseForTag(pyxb.namespace.ExpandedName(None, 'Swst_changes_values')), pyxb.utils.utility.Location('aresysTypes.xsd', 935, 10))
    st_9 = fac.State(symbol, is_initial=False, final_update=final_update, is_unordered_catenation=False)
    states.append(st_9)
    final_update = None
    symbol = pyxb.binding.content.ElementUse(AcquisitionTimelineType._UseForTag(pyxb.namespace.ExpandedName(None, 'noise_packets_number')), pyxb.utils.utility.Location('aresysTypes.xsd', 954, 10))
    st_10 = fac.State(symbol, is_initial=False, final_update=final_update, is_unordered_catenation=False)
    states.append(st_10)
    final_update = None
    symbol = pyxb.binding.content.ElementUse(AcquisitionTimelineType._UseForTag(pyxb.namespace.ExpandedName(None, 'noise_packets_azimuthtimes')), pyxb.utils.utility.Location('aresysTypes.xsd', 959, 10))
    st_11 = fac.State(symbol, is_initial=False, final_update=final_update, is_unordered_catenation=False)
    states.append(st_11)
    final_update = None
    symbol = pyxb.binding.content.ElementUse(AcquisitionTimelineType._UseForTag(pyxb.namespace.ExpandedName(None, 'Internal_calibration_number')), pyxb.utils.utility.Location('aresysTypes.xsd', 978, 10))
    st_12 = fac.State(symbol, is_initial=False, final_update=final_update, is_unordered_catenation=False)
    states.append(st_12)
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(AcquisitionTimelineType._UseForTag(pyxb.namespace.ExpandedName(None, 'Internal_calibration_azimuthtimes')), pyxb.utils.utility.Location('aresysTypes.xsd', 983, 10))
    st_13 = fac.State(symbol, is_initial=False, final_update=final_update, is_unordered_catenation=False)
    states.append(st_13)
    final_update = None
    symbol = pyxb.binding.content.ElementUse(AcquisitionTimelineType._UseForTag(pyxb.namespace.ExpandedName(None, 'Swl_changes_number')), pyxb.utils.utility.Location('aresysTypes.xsd', 1003, 12))
    st_14 = fac.State(symbol, is_initial=False, final_update=final_update, is_unordered_catenation=False)
    states.append(st_14)
    final_update = None
    symbol = pyxb.binding.content.ElementUse(AcquisitionTimelineType._UseForTag(pyxb.namespace.ExpandedName(None, 'Swl_changes_azimuthtimes')), pyxb.utils.utility.Location('aresysTypes.xsd', 1008, 12))
    st_15 = fac.State(symbol, is_initial=False, final_update=final_update, is_unordered_catenation=False)
    states.append(st_15)
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_5, False))
    symbol = pyxb.binding.content.ElementUse(AcquisitionTimelineType._UseForTag(pyxb.namespace.ExpandedName(None, 'Swl_changes_values')), pyxb.utils.utility.Location('aresysTypes.xsd', 1027, 12))
    st_16 = fac.State(symbol, is_initial=False, final_update=final_update, is_unordered_catenation=False)
    states.append(st_16)
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_6, False))
    symbol = pyxb.binding.content.ElementUse(AcquisitionTimelineType._UseForTag(pyxb.namespace.ExpandedName(None, 'ChirpPeriod')), pyxb.utils.utility.Location('aresysTypes.xsd', 1047, 10))
    st_17 = fac.State(symbol, is_initial=False, final_update=final_update, is_unordered_catenation=False)
    states.append(st_17)
    transitions = []
    transitions.append(fac.Transition(st_1, [
         ]))
    st_0._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_2, [
         ]))
    transitions.append(fac.Transition(st_3, [
         ]))
    transitions.append(fac.Transition(st_4, [
         ]))
    transitions.append(fac.Transition(st_5, [
         ]))
    transitions.append(fac.Transition(st_6, [
         ]))
    transitions.append(fac.Transition(st_7, [
         ]))
    st_1._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_2, [
        fac.UpdateInstruction(cc_0, True) ]))
    transitions.append(fac.Transition(st_3, [
        fac.UpdateInstruction(cc_0, False) ]))
    transitions.append(fac.Transition(st_4, [
        fac.UpdateInstruction(cc_0, False) ]))
    transitions.append(fac.Transition(st_5, [
        fac.UpdateInstruction(cc_0, False) ]))
    transitions.append(fac.Transition(st_6, [
        fac.UpdateInstruction(cc_0, False) ]))
    transitions.append(fac.Transition(st_7, [
        fac.UpdateInstruction(cc_0, False) ]))
    st_2._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_3, [
        fac.UpdateInstruction(cc_1, True) ]))
    transitions.append(fac.Transition(st_4, [
        fac.UpdateInstruction(cc_1, False) ]))
    transitions.append(fac.Transition(st_5, [
        fac.UpdateInstruction(cc_1, False) ]))
    transitions.append(fac.Transition(st_6, [
        fac.UpdateInstruction(cc_1, False) ]))
    transitions.append(fac.Transition(st_7, [
        fac.UpdateInstruction(cc_1, False) ]))
    st_3._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_4, [
        fac.UpdateInstruction(cc_2, True) ]))
    transitions.append(fac.Transition(st_5, [
        fac.UpdateInstruction(cc_2, False) ]))
    transitions.append(fac.Transition(st_6, [
        fac.UpdateInstruction(cc_2, False) ]))
    transitions.append(fac.Transition(st_7, [
        fac.UpdateInstruction(cc_2, False) ]))
    st_4._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_5, [
        fac.UpdateInstruction(cc_3, True) ]))
    transitions.append(fac.Transition(st_6, [
        fac.UpdateInstruction(cc_3, False) ]))
    transitions.append(fac.Transition(st_7, [
        fac.UpdateInstruction(cc_3, False) ]))
    st_5._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_6, [
        fac.UpdateInstruction(cc_4, True) ]))
    transitions.append(fac.Transition(st_7, [
        fac.UpdateInstruction(cc_4, False) ]))
    st_6._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_8, [
         ]))
    st_7._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_9, [
         ]))
    st_8._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_10, [
         ]))
    st_9._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_11, [
         ]))
    st_10._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_12, [
         ]))
    st_11._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_13, [
         ]))
    st_12._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_14, [
         ]))
    transitions.append(fac.Transition(st_17, [
         ]))
    st_13._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_15, [
         ]))
    st_14._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_16, [
         ]))
    st_15._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_14, [
        fac.UpdateInstruction(cc_5, True) ]))
    transitions.append(fac.Transition(st_17, [
        fac.UpdateInstruction(cc_5, False) ]))
    st_16._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_17, [
        fac.UpdateInstruction(cc_6, True) ]))
    st_17._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=None)
AcquisitionTimelineType._Automaton = _BuildAutomaton_98()




RasterInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'FileName'), pyxb.binding.datatypes.string, scope=RasterInfoType, documentation='Name of the associated binary file', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1063, 10)))

RasterInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'Lines'), pyxb.binding.datatypes.unsignedInt, scope=RasterInfoType, documentation='Total number of lines (azimuth) of the image', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1068, 10)))

RasterInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'Samples'), pyxb.binding.datatypes.unsignedInt, scope=RasterInfoType, documentation='Total number of samples (range) of the image', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1073, 10)))

RasterInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'HeaderOffsetBytes'), pyxb.binding.datatypes.unsignedInt, scope=RasterInfoType, documentation='Number of bytes at the beginning of the file containing the header information', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1078, 10)))

RasterInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'RowPrefixBytes'), pyxb.binding.datatypes.unsignedInt, scope=RasterInfoType, documentation='Number of bytes at the beginning of each line containing header information', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1083, 10)))

RasterInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'ByteOrder'), Endianity, scope=RasterInfoType, documentation='Endianity: BIGENDIAN or LITTLEENDIAN', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1088, 10)))

RasterInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'CellType'), CellTypeVerboseType, scope=RasterInfoType, documentation='Byte format type: FLOAT_COMPLEX, FLOAT32, DOUBLE_COMPLEX, FLOAT64, INT16, SHORT_COMPLEX, INT32, INT_COMPLEX, INT8, INT8_COMPLEX', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1093, 10)))

RasterInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'LinesStep'), doubleWithUnit, scope=RasterInfoType, documentation='Azimuth sampling step [s]', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1098, 10)))

RasterInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'SamplesStep'), doubleWithUnit, scope=RasterInfoType, documentation='Range sampling step [s]', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1103, 10)))

RasterInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'LinesStart'), stringWithUnit, scope=RasterInfoType, documentation='Azimuth absolute start time [Utc]', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1108, 10)))

RasterInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'SamplesStart'), stringWithUnit, scope=RasterInfoType, documentation='Range absolute start time [s]', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1113, 10)))

RasterInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'RasterFormat'), RasterFormatType, scope=RasterInfoType, documentation='Raster Format of the Data ( Default value is ARESYS_RASTER )', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1118, 10)))

RasterInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'InvalidValue'), DCOMPLEX, scope=RasterInfoType, location=pyxb.utils.utility.Location('aresysTypes.xsd', 1123, 10)))

def _BuildAutomaton_100 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_100
    del _BuildAutomaton_100
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(RasterInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'FileName')), pyxb.utils.utility.Location('aresysTypes.xsd', 1063, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_101 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_101
    del _BuildAutomaton_101
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(RasterInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'Lines')), pyxb.utils.utility.Location('aresysTypes.xsd', 1068, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_102 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_102
    del _BuildAutomaton_102
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(RasterInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'Samples')), pyxb.utils.utility.Location('aresysTypes.xsd', 1073, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_103 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_103
    del _BuildAutomaton_103
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(RasterInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'HeaderOffsetBytes')), pyxb.utils.utility.Location('aresysTypes.xsd', 1078, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_104 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_104
    del _BuildAutomaton_104
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(RasterInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'RowPrefixBytes')), pyxb.utils.utility.Location('aresysTypes.xsd', 1083, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_105 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_105
    del _BuildAutomaton_105
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(RasterInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'ByteOrder')), pyxb.utils.utility.Location('aresysTypes.xsd', 1088, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_106 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_106
    del _BuildAutomaton_106
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(RasterInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'CellType')), pyxb.utils.utility.Location('aresysTypes.xsd', 1093, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_107 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_107
    del _BuildAutomaton_107
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(RasterInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'LinesStep')), pyxb.utils.utility.Location('aresysTypes.xsd', 1098, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_108 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_108
    del _BuildAutomaton_108
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(RasterInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'SamplesStep')), pyxb.utils.utility.Location('aresysTypes.xsd', 1103, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_109 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_109
    del _BuildAutomaton_109
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(RasterInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'LinesStart')), pyxb.utils.utility.Location('aresysTypes.xsd', 1108, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_110 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_110
    del _BuildAutomaton_110
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(RasterInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'SamplesStart')), pyxb.utils.utility.Location('aresysTypes.xsd', 1113, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_111 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_111
    del _BuildAutomaton_111
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 1118, 10))
    counters.add(cc_0)
    states = []
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_0, False))
    symbol = pyxb.binding.content.ElementUse(RasterInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'RasterFormat')), pyxb.utils.utility.Location('aresysTypes.xsd', 1118, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    transitions.append(fac.Transition(st_0, [
        fac.UpdateInstruction(cc_0, True) ]))
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, True, containing_state=st_0)

def _BuildAutomaton_112 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_112
    del _BuildAutomaton_112
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 1123, 10))
    counters.add(cc_0)
    states = []
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_0, False))
    symbol = pyxb.binding.content.ElementUse(RasterInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'InvalidValue')), pyxb.utils.utility.Location('aresysTypes.xsd', 1123, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    transitions.append(fac.Transition(st_0, [
        fac.UpdateInstruction(cc_0, True) ]))
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, True, containing_state=st_0)

def _BuildAutomaton_99 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_99
    del _BuildAutomaton_99
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 1118, 10))
    counters.add(cc_0)
    cc_1 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 1123, 10))
    counters.add(cc_1)
    states = []
    sub_automata = []
    sub_automata.append(_BuildAutomaton_100())
    sub_automata.append(_BuildAutomaton_101())
    sub_automata.append(_BuildAutomaton_102())
    sub_automata.append(_BuildAutomaton_103())
    sub_automata.append(_BuildAutomaton_104())
    sub_automata.append(_BuildAutomaton_105())
    sub_automata.append(_BuildAutomaton_106())
    sub_automata.append(_BuildAutomaton_107())
    sub_automata.append(_BuildAutomaton_108())
    sub_automata.append(_BuildAutomaton_109())
    sub_automata.append(_BuildAutomaton_110())
    sub_automata.append(_BuildAutomaton_111())
    sub_automata.append(_BuildAutomaton_112())
    final_update = set()
    symbol = pyxb.utils.utility.Location('aresysTypes.xsd', 1062, 8)
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=True)
    st_0._set_subAutomata(*sub_automata)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=None)
RasterInfoType._Automaton = _BuildAutomaton_99()




SamplingConstantsType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'frg_hz'), CTD_ANON_72, scope=SamplingConstantsType, documentation='Range sampling frequency [Hz]', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1135, 10)))

SamplingConstantsType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'Brg_hz'), CTD_ANON_73, scope=SamplingConstantsType, documentation='Range bandwidth [Hz]', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1147, 10)))

SamplingConstantsType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'faz_hz'), CTD_ANON_74, scope=SamplingConstantsType, documentation='Azimuth sampling frequency [Hz]', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1159, 10)))

SamplingConstantsType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'Baz_hz'), CTD_ANON_75, scope=SamplingConstantsType, documentation='Azimuth bandwidth [Hz]', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1171, 10)))

def _BuildAutomaton_114 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_114
    del _BuildAutomaton_114
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(SamplingConstantsType._UseForTag(pyxb.namespace.ExpandedName(None, 'frg_hz')), pyxb.utils.utility.Location('aresysTypes.xsd', 1135, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_115 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_115
    del _BuildAutomaton_115
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(SamplingConstantsType._UseForTag(pyxb.namespace.ExpandedName(None, 'Brg_hz')), pyxb.utils.utility.Location('aresysTypes.xsd', 1147, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_116 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_116
    del _BuildAutomaton_116
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(SamplingConstantsType._UseForTag(pyxb.namespace.ExpandedName(None, 'faz_hz')), pyxb.utils.utility.Location('aresysTypes.xsd', 1159, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_117 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_117
    del _BuildAutomaton_117
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(SamplingConstantsType._UseForTag(pyxb.namespace.ExpandedName(None, 'Baz_hz')), pyxb.utils.utility.Location('aresysTypes.xsd', 1171, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_113 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_113
    del _BuildAutomaton_113
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    sub_automata = []
    sub_automata.append(_BuildAutomaton_114())
    sub_automata.append(_BuildAutomaton_115())
    sub_automata.append(_BuildAutomaton_116())
    sub_automata.append(_BuildAutomaton_117())
    final_update = set()
    symbol = pyxb.utils.utility.Location('aresysTypes.xsd', 1134, 8)
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=True)
    st_0._set_subAutomata(*sub_automata)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=None)
SamplingConstantsType._Automaton = _BuildAutomaton_113()




DataStatisticsType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'NumSamples'), CTD_ANON_76, scope=DataStatisticsType, documentation='Number of samples analyzed', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1194, 10)))

DataStatisticsType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'MaxI'), CTD_ANON_77, scope=DataStatisticsType, documentation='Max of I (real) samples', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1206, 10)))

DataStatisticsType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'MinI'), CTD_ANON_78, scope=DataStatisticsType, documentation='Min of I (real) samples', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1218, 10)))

DataStatisticsType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'MaxQ'), CTD_ANON_79, scope=DataStatisticsType, documentation='Max of Q (imaginary) samples', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1230, 10)))

DataStatisticsType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'MinQ'), CTD_ANON_80, scope=DataStatisticsType, documentation='Min of Q (imaginary) samples', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1242, 10)))

DataStatisticsType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'SumI'), CTD_ANON_81, scope=DataStatisticsType, documentation='Sum of I (real) samples', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1254, 10)))

DataStatisticsType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'SumQ'), CTD_ANON_82, scope=DataStatisticsType, documentation='Sum of Q (imaginary) samples', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1266, 10)))

DataStatisticsType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'Sum2I'), CTD_ANON_83, scope=DataStatisticsType, documentation='Square Sum of I (real) samples', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1278, 10)))

DataStatisticsType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'Sum2Q'), CTD_ANON_84, scope=DataStatisticsType, documentation='Square Sum of Q (imaginary) samples', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1290, 10)))

DataStatisticsType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'StdDevI'), CTD_ANON_85, scope=DataStatisticsType, documentation='Standard Deviation of I (real) samples', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1302, 10)))

DataStatisticsType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'StdDevQ'), CTD_ANON_86, scope=DataStatisticsType, documentation='Standard Deviation of Q (imaginary) samples', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1314, 10)))

DataStatisticsType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'StatisticsList'), CTD_ANON_24, scope=DataStatisticsType, location=pyxb.utils.utility.Location('aresysTypes.xsd', 1326, 10)))

def _BuildAutomaton_119 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_119
    del _BuildAutomaton_119
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(DataStatisticsType._UseForTag(pyxb.namespace.ExpandedName(None, 'NumSamples')), pyxb.utils.utility.Location('aresysTypes.xsd', 1194, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_120 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_120
    del _BuildAutomaton_120
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(DataStatisticsType._UseForTag(pyxb.namespace.ExpandedName(None, 'MaxI')), pyxb.utils.utility.Location('aresysTypes.xsd', 1206, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_121 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_121
    del _BuildAutomaton_121
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(DataStatisticsType._UseForTag(pyxb.namespace.ExpandedName(None, 'MinI')), pyxb.utils.utility.Location('aresysTypes.xsd', 1218, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_122 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_122
    del _BuildAutomaton_122
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(DataStatisticsType._UseForTag(pyxb.namespace.ExpandedName(None, 'MaxQ')), pyxb.utils.utility.Location('aresysTypes.xsd', 1230, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_123 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_123
    del _BuildAutomaton_123
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(DataStatisticsType._UseForTag(pyxb.namespace.ExpandedName(None, 'MinQ')), pyxb.utils.utility.Location('aresysTypes.xsd', 1242, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_124 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_124
    del _BuildAutomaton_124
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(DataStatisticsType._UseForTag(pyxb.namespace.ExpandedName(None, 'SumI')), pyxb.utils.utility.Location('aresysTypes.xsd', 1254, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_125 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_125
    del _BuildAutomaton_125
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(DataStatisticsType._UseForTag(pyxb.namespace.ExpandedName(None, 'SumQ')), pyxb.utils.utility.Location('aresysTypes.xsd', 1266, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_126 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_126
    del _BuildAutomaton_126
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(DataStatisticsType._UseForTag(pyxb.namespace.ExpandedName(None, 'Sum2I')), pyxb.utils.utility.Location('aresysTypes.xsd', 1278, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_127 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_127
    del _BuildAutomaton_127
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(DataStatisticsType._UseForTag(pyxb.namespace.ExpandedName(None, 'Sum2Q')), pyxb.utils.utility.Location('aresysTypes.xsd', 1290, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_128 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_128
    del _BuildAutomaton_128
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(DataStatisticsType._UseForTag(pyxb.namespace.ExpandedName(None, 'StdDevI')), pyxb.utils.utility.Location('aresysTypes.xsd', 1302, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_129 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_129
    del _BuildAutomaton_129
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(DataStatisticsType._UseForTag(pyxb.namespace.ExpandedName(None, 'StdDevQ')), pyxb.utils.utility.Location('aresysTypes.xsd', 1314, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_130 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_130
    del _BuildAutomaton_130
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 1326, 10))
    counters.add(cc_0)
    states = []
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_0, False))
    symbol = pyxb.binding.content.ElementUse(DataStatisticsType._UseForTag(pyxb.namespace.ExpandedName(None, 'StatisticsList')), pyxb.utils.utility.Location('aresysTypes.xsd', 1326, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    transitions.append(fac.Transition(st_0, [
        fac.UpdateInstruction(cc_0, True) ]))
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, True, containing_state=st_0)

def _BuildAutomaton_118 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_118
    del _BuildAutomaton_118
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 1326, 10))
    counters.add(cc_0)
    states = []
    sub_automata = []
    sub_automata.append(_BuildAutomaton_119())
    sub_automata.append(_BuildAutomaton_120())
    sub_automata.append(_BuildAutomaton_121())
    sub_automata.append(_BuildAutomaton_122())
    sub_automata.append(_BuildAutomaton_123())
    sub_automata.append(_BuildAutomaton_124())
    sub_automata.append(_BuildAutomaton_125())
    sub_automata.append(_BuildAutomaton_126())
    sub_automata.append(_BuildAutomaton_127())
    sub_automata.append(_BuildAutomaton_128())
    sub_automata.append(_BuildAutomaton_129())
    sub_automata.append(_BuildAutomaton_130())
    final_update = set()
    symbol = pyxb.utils.utility.Location('aresysTypes.xsd', 1193, 8)
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=True)
    st_0._set_subAutomata(*sub_automata)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=None)
DataStatisticsType._Automaton = _BuildAutomaton_118()




DataBlockStatisticsType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'NumSamples'), CTD_ANON_87, scope=DataBlockStatisticsType, documentation='Number of samples analyzed', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1344, 10)))

DataBlockStatisticsType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'MaxI'), CTD_ANON_88, scope=DataBlockStatisticsType, documentation='Max of I (real) samples', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1356, 10)))

DataBlockStatisticsType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'MinI'), CTD_ANON_89, scope=DataBlockStatisticsType, documentation='Min of I (real) samples', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1368, 10)))

DataBlockStatisticsType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'MaxQ'), CTD_ANON_90, scope=DataBlockStatisticsType, documentation='Max of Q (imaginary) samples', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1380, 10)))

DataBlockStatisticsType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'MinQ'), CTD_ANON_91, scope=DataBlockStatisticsType, documentation='Min of Q (imaginary) samples', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1392, 10)))

DataBlockStatisticsType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'SumI'), CTD_ANON_92, scope=DataBlockStatisticsType, documentation='Sum of I (real) samples', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1404, 10)))

DataBlockStatisticsType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'SumQ'), CTD_ANON_93, scope=DataBlockStatisticsType, documentation='Sum of Q (imaginary) samples', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1416, 10)))

DataBlockStatisticsType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'Sum2I'), CTD_ANON_94, scope=DataBlockStatisticsType, documentation='Square Sum of I (real) samples', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1428, 10)))

DataBlockStatisticsType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'Sum2Q'), CTD_ANON_95, scope=DataBlockStatisticsType, documentation='Square Sum of Q (imaginary) samples', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1440, 10)))

def _BuildAutomaton_132 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_132
    del _BuildAutomaton_132
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(DataBlockStatisticsType._UseForTag(pyxb.namespace.ExpandedName(None, 'NumSamples')), pyxb.utils.utility.Location('aresysTypes.xsd', 1344, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_133 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_133
    del _BuildAutomaton_133
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(DataBlockStatisticsType._UseForTag(pyxb.namespace.ExpandedName(None, 'MaxI')), pyxb.utils.utility.Location('aresysTypes.xsd', 1356, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_134 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_134
    del _BuildAutomaton_134
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(DataBlockStatisticsType._UseForTag(pyxb.namespace.ExpandedName(None, 'MinI')), pyxb.utils.utility.Location('aresysTypes.xsd', 1368, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_135 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_135
    del _BuildAutomaton_135
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(DataBlockStatisticsType._UseForTag(pyxb.namespace.ExpandedName(None, 'MaxQ')), pyxb.utils.utility.Location('aresysTypes.xsd', 1380, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_136 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_136
    del _BuildAutomaton_136
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(DataBlockStatisticsType._UseForTag(pyxb.namespace.ExpandedName(None, 'MinQ')), pyxb.utils.utility.Location('aresysTypes.xsd', 1392, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_137 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_137
    del _BuildAutomaton_137
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(DataBlockStatisticsType._UseForTag(pyxb.namespace.ExpandedName(None, 'SumI')), pyxb.utils.utility.Location('aresysTypes.xsd', 1404, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_138 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_138
    del _BuildAutomaton_138
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(DataBlockStatisticsType._UseForTag(pyxb.namespace.ExpandedName(None, 'SumQ')), pyxb.utils.utility.Location('aresysTypes.xsd', 1416, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_139 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_139
    del _BuildAutomaton_139
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(DataBlockStatisticsType._UseForTag(pyxb.namespace.ExpandedName(None, 'Sum2I')), pyxb.utils.utility.Location('aresysTypes.xsd', 1428, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_140 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_140
    del _BuildAutomaton_140
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(DataBlockStatisticsType._UseForTag(pyxb.namespace.ExpandedName(None, 'Sum2Q')), pyxb.utils.utility.Location('aresysTypes.xsd', 1440, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_131 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_131
    del _BuildAutomaton_131
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    sub_automata = []
    sub_automata.append(_BuildAutomaton_132())
    sub_automata.append(_BuildAutomaton_133())
    sub_automata.append(_BuildAutomaton_134())
    sub_automata.append(_BuildAutomaton_135())
    sub_automata.append(_BuildAutomaton_136())
    sub_automata.append(_BuildAutomaton_137())
    sub_automata.append(_BuildAutomaton_138())
    sub_automata.append(_BuildAutomaton_139())
    sub_automata.append(_BuildAutomaton_140())
    final_update = set()
    symbol = pyxb.utils.utility.Location('aresysTypes.xsd', 1343, 8)
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=True)
    st_0._set_subAutomata(*sub_automata)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=None)
DataBlockStatisticsType._Automaton = _BuildAutomaton_131()




polyType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'pol'), CTD_ANON_25, scope=polyType, documentation='Polynomial coefficients: const, rg, az, az*rg, rg^2, rg^3, rg^4 [Optional: rg^5 rg^6 .... rg^N]', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1465, 10)))

polyType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'trg0_s'), CTD_ANON_97, scope=polyType, documentation='Polynomial range reference time [s]', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1484, 10)))

polyType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'taz0_Utc'), CTD_ANON_98, scope=polyType, documentation='Polynomial azimuth reference time [Utc]', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1496, 10)))

def _BuildAutomaton_141 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_141
    del _BuildAutomaton_141
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = None
    symbol = pyxb.binding.content.ElementUse(polyType._UseForTag(pyxb.namespace.ExpandedName(None, 'pol')), pyxb.utils.utility.Location('aresysTypes.xsd', 1465, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    final_update = None
    symbol = pyxb.binding.content.ElementUse(polyType._UseForTag(pyxb.namespace.ExpandedName(None, 'trg0_s')), pyxb.utils.utility.Location('aresysTypes.xsd', 1484, 10))
    st_1 = fac.State(symbol, is_initial=False, final_update=final_update, is_unordered_catenation=False)
    states.append(st_1)
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(polyType._UseForTag(pyxb.namespace.ExpandedName(None, 'taz0_Utc')), pyxb.utils.utility.Location('aresysTypes.xsd', 1496, 10))
    st_2 = fac.State(symbol, is_initial=False, final_update=final_update, is_unordered_catenation=False)
    states.append(st_2)
    transitions = []
    transitions.append(fac.Transition(st_1, [
         ]))
    st_0._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_2, [
         ]))
    st_1._set_transitionSet(transitions)
    transitions = []
    st_2._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=None)
polyType._Automaton = _BuildAutomaton_141()




polyCoregType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'polRg'), CTD_ANON_26, scope=polyCoregType, documentation='Polynomial coefficients of Range', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1519, 10)))

polyCoregType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'polAz'), CTD_ANON_27, scope=polyCoregType, documentation='Polynomial coefficients of Azimuth', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1538, 10)))

polyCoregType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'trg0_s'), CTD_ANON_101, scope=polyCoregType, documentation='Polynomial range reference time [s]', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1557, 10)))

polyCoregType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'taz0_Utc'), CTD_ANON_102, scope=polyCoregType, documentation='Polynomial azimuth reference time [Utc]', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1569, 10)))

def _BuildAutomaton_142 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_142
    del _BuildAutomaton_142
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = None
    symbol = pyxb.binding.content.ElementUse(polyCoregType._UseForTag(pyxb.namespace.ExpandedName(None, 'polRg')), pyxb.utils.utility.Location('aresysTypes.xsd', 1519, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    final_update = None
    symbol = pyxb.binding.content.ElementUse(polyCoregType._UseForTag(pyxb.namespace.ExpandedName(None, 'polAz')), pyxb.utils.utility.Location('aresysTypes.xsd', 1538, 10))
    st_1 = fac.State(symbol, is_initial=False, final_update=final_update, is_unordered_catenation=False)
    states.append(st_1)
    final_update = None
    symbol = pyxb.binding.content.ElementUse(polyCoregType._UseForTag(pyxb.namespace.ExpandedName(None, 'trg0_s')), pyxb.utils.utility.Location('aresysTypes.xsd', 1557, 10))
    st_2 = fac.State(symbol, is_initial=False, final_update=final_update, is_unordered_catenation=False)
    states.append(st_2)
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(polyCoregType._UseForTag(pyxb.namespace.ExpandedName(None, 'taz0_Utc')), pyxb.utils.utility.Location('aresysTypes.xsd', 1569, 10))
    st_3 = fac.State(symbol, is_initial=False, final_update=final_update, is_unordered_catenation=False)
    states.append(st_3)
    transitions = []
    transitions.append(fac.Transition(st_1, [
         ]))
    st_0._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_2, [
         ]))
    st_1._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_3, [
         ]))
    st_2._set_transitionSet(transitions)
    transitions = []
    st_3._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=None)
polyCoregType._Automaton = _BuildAutomaton_142()




GroundCornersPointsType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'EastingGridSize'), CTD_ANON_103, scope=GroundCornersPointsType, location=pyxb.utils.utility.Location('aresysTypes.xsd', 1589, 10)))

GroundCornersPointsType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'NorthingGridSize'), CTD_ANON_104, scope=GroundCornersPointsType, location=pyxb.utils.utility.Location('aresysTypes.xsd', 1598, 10)))

GroundCornersPointsType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'NorthWest'), CTD_ANON_28, scope=GroundCornersPointsType, location=pyxb.utils.utility.Location('aresysTypes.xsd', 1607, 10)))

GroundCornersPointsType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'NorthEast'), CTD_ANON_29, scope=GroundCornersPointsType, location=pyxb.utils.utility.Location('aresysTypes.xsd', 1614, 10)))

GroundCornersPointsType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'SouthWest'), CTD_ANON_30, scope=GroundCornersPointsType, location=pyxb.utils.utility.Location('aresysTypes.xsd', 1621, 10)))

GroundCornersPointsType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'SouthEast'), CTD_ANON_31, scope=GroundCornersPointsType, location=pyxb.utils.utility.Location('aresysTypes.xsd', 1628, 10)))

GroundCornersPointsType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'Center'), CTD_ANON_32, scope=GroundCornersPointsType, location=pyxb.utils.utility.Location('aresysTypes.xsd', 1635, 10)))

def _BuildAutomaton_144 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_144
    del _BuildAutomaton_144
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(GroundCornersPointsType._UseForTag(pyxb.namespace.ExpandedName(None, 'EastingGridSize')), pyxb.utils.utility.Location('aresysTypes.xsd', 1589, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_145 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_145
    del _BuildAutomaton_145
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(GroundCornersPointsType._UseForTag(pyxb.namespace.ExpandedName(None, 'NorthingGridSize')), pyxb.utils.utility.Location('aresysTypes.xsd', 1598, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_146 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_146
    del _BuildAutomaton_146
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(GroundCornersPointsType._UseForTag(pyxb.namespace.ExpandedName(None, 'NorthWest')), pyxb.utils.utility.Location('aresysTypes.xsd', 1607, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_147 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_147
    del _BuildAutomaton_147
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(GroundCornersPointsType._UseForTag(pyxb.namespace.ExpandedName(None, 'NorthEast')), pyxb.utils.utility.Location('aresysTypes.xsd', 1614, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_148 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_148
    del _BuildAutomaton_148
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(GroundCornersPointsType._UseForTag(pyxb.namespace.ExpandedName(None, 'SouthWest')), pyxb.utils.utility.Location('aresysTypes.xsd', 1621, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_149 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_149
    del _BuildAutomaton_149
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(GroundCornersPointsType._UseForTag(pyxb.namespace.ExpandedName(None, 'SouthEast')), pyxb.utils.utility.Location('aresysTypes.xsd', 1628, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_150 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_150
    del _BuildAutomaton_150
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(GroundCornersPointsType._UseForTag(pyxb.namespace.ExpandedName(None, 'Center')), pyxb.utils.utility.Location('aresysTypes.xsd', 1635, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_143 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_143
    del _BuildAutomaton_143
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    sub_automata = []
    sub_automata.append(_BuildAutomaton_144())
    sub_automata.append(_BuildAutomaton_145())
    sub_automata.append(_BuildAutomaton_146())
    sub_automata.append(_BuildAutomaton_147())
    sub_automata.append(_BuildAutomaton_148())
    sub_automata.append(_BuildAutomaton_149())
    sub_automata.append(_BuildAutomaton_150())
    final_update = set()
    symbol = pyxb.utils.utility.Location('aresysTypes.xsd', 1588, 8)
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=True)
    st_0._set_subAutomata(*sub_automata)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=None)
GroundCornersPointsType._Automaton = _BuildAutomaton_143()




BurstInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'NumberOfBursts'), pyxb.binding.datatypes.unsignedInt, scope=BurstInfoType, documentation='Number of bursts in the swath', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1666, 10)))

BurstInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'LinesPerBurst'), pyxb.binding.datatypes.unsignedInt, scope=BurstInfoType, documentation='Number of lines in each burst', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1672, 12)))

BurstInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'LinesPerBurstChangeList'), CTD_ANON_33, scope=BurstInfoType, location=pyxb.utils.utility.Location('aresysTypes.xsd', 1677, 12)))

BurstInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'BurstRepetitionFrequency'), doubleWithUnit, scope=BurstInfoType, documentation='Burst repetition frequency [Hz]', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1693, 10)))

BurstInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'Burst'), BurstType, scope=BurstInfoType, documentation='Time coordinates of each burst', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1698, 10)))

def _BuildAutomaton_151 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_151
    del _BuildAutomaton_151
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = None
    symbol = pyxb.binding.content.ElementUse(BurstInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'NumberOfBursts')), pyxb.utils.utility.Location('aresysTypes.xsd', 1666, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    final_update = None
    symbol = pyxb.binding.content.ElementUse(BurstInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'LinesPerBurst')), pyxb.utils.utility.Location('aresysTypes.xsd', 1672, 12))
    st_1 = fac.State(symbol, is_initial=False, final_update=final_update, is_unordered_catenation=False)
    states.append(st_1)
    final_update = None
    symbol = pyxb.binding.content.ElementUse(BurstInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'LinesPerBurstChangeList')), pyxb.utils.utility.Location('aresysTypes.xsd', 1677, 12))
    st_2 = fac.State(symbol, is_initial=False, final_update=final_update, is_unordered_catenation=False)
    states.append(st_2)
    final_update = None
    symbol = pyxb.binding.content.ElementUse(BurstInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'BurstRepetitionFrequency')), pyxb.utils.utility.Location('aresysTypes.xsd', 1693, 10))
    st_3 = fac.State(symbol, is_initial=False, final_update=final_update, is_unordered_catenation=False)
    states.append(st_3)
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(BurstInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'Burst')), pyxb.utils.utility.Location('aresysTypes.xsd', 1698, 10))
    st_4 = fac.State(symbol, is_initial=False, final_update=final_update, is_unordered_catenation=False)
    states.append(st_4)
    transitions = []
    transitions.append(fac.Transition(st_1, [
         ]))
    transitions.append(fac.Transition(st_2, [
         ]))
    st_0._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_3, [
         ]))
    st_1._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_3, [
         ]))
    st_2._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_4, [
         ]))
    st_3._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_4, [
         ]))
    st_4._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=None)
BurstInfoType._Automaton = _BuildAutomaton_151()




AntennaInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'SensorName'), SensorNamesType, scope=AntennaInfoType, documentation='Sensor name: ASAR, PALSAR, ...', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1714, 10)))

AntennaInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'AcquisitionMode'), AcquisitionModeType, scope=AntennaInfoType, documentation='Acquisition mode: STRIPMAP, TOPSAR, ...', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1719, 10)))

AntennaInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'BeamName'), pyxb.binding.datatypes.string, scope=AntennaInfoType, documentation='Acquisition beam name', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1724, 10)))

AntennaInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'Polarization'), PolarizationType, scope=AntennaInfoType, documentation='Antenna polarization (H/H, H/V...)', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1729, 10)))

AntennaInfoType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'LinesPerPattern'), pyxb.binding.datatypes.unsignedInt, scope=AntennaInfoType, documentation='Contains number of lines for each pattern', location=pyxb.utils.utility.Location('aresysTypes.xsd', 1734, 10)))

def _BuildAutomaton_152 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_152
    del _BuildAutomaton_152
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresysTypes.xsd', 1734, 10))
    counters.add(cc_0)
    states = []
    final_update = None
    symbol = pyxb.binding.content.ElementUse(AntennaInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'SensorName')), pyxb.utils.utility.Location('aresysTypes.xsd', 1714, 10))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    final_update = None
    symbol = pyxb.binding.content.ElementUse(AntennaInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'AcquisitionMode')), pyxb.utils.utility.Location('aresysTypes.xsd', 1719, 10))
    st_1 = fac.State(symbol, is_initial=False, final_update=final_update, is_unordered_catenation=False)
    states.append(st_1)
    final_update = None
    symbol = pyxb.binding.content.ElementUse(AntennaInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'BeamName')), pyxb.utils.utility.Location('aresysTypes.xsd', 1724, 10))
    st_2 = fac.State(symbol, is_initial=False, final_update=final_update, is_unordered_catenation=False)
    states.append(st_2)
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(AntennaInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'Polarization')), pyxb.utils.utility.Location('aresysTypes.xsd', 1729, 10))
    st_3 = fac.State(symbol, is_initial=False, final_update=final_update, is_unordered_catenation=False)
    states.append(st_3)
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_0, False))
    symbol = pyxb.binding.content.ElementUse(AntennaInfoType._UseForTag(pyxb.namespace.ExpandedName(None, 'LinesPerPattern')), pyxb.utils.utility.Location('aresysTypes.xsd', 1734, 10))
    st_4 = fac.State(symbol, is_initial=False, final_update=final_update, is_unordered_catenation=False)
    states.append(st_4)
    transitions = []
    transitions.append(fac.Transition(st_1, [
         ]))
    st_0._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_2, [
         ]))
    st_1._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_3, [
         ]))
    st_2._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_4, [
         ]))
    st_3._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_4, [
        fac.UpdateInstruction(cc_0, True) ]))
    st_4._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=None)
AntennaInfoType._Automaton = _BuildAutomaton_152()

