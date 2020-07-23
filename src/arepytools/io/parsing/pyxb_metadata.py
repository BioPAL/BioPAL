# arepytools/io/parsing/pyxb_metadata.py
# -*- coding: utf-8 -*-
# PyXB bindings for NM:e92452c8d3e28a9e27abfc9994d2007779e7f4c9
# Generated 2021-02-17 11:26:26.605622 by PyXB version 1.2.6 using Python 3.8.5.final.0
# Namespace AbsentNamespace0

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
from . import pyxb_metadata_types as _ImportedBinding_pyxb_metadata_types

# NOTE: All namespace declarations are reserved within the binding
Namespace = pyxb.namespace.CreateAbsentNamespace()
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


# Complex type AresysXmlDocType with content type ELEMENT_ONLY
class AresysXmlDocType (pyxb.binding.basis.complexTypeDefinition):
    """Complex type AresysXmlDocType with content type ELEMENT_ONLY"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'AresysXmlDocType')
    _XSDLocation = pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 5, 2)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Element NumberOfChannels uses Python identifier NumberOfChannels
    __NumberOfChannels = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'NumberOfChannels'), 'NumberOfChannels', '__AbsentNamespace0_AresysXmlDocType_NumberOfChannels', False, pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 7, 6), )

    
    NumberOfChannels = property(__NumberOfChannels.value, __NumberOfChannels.set, None, None)

    
    # Element VersionNumber uses Python identifier VersionNumber
    __VersionNumber = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'VersionNumber'), 'VersionNumber', '__AbsentNamespace0_AresysXmlDocType_VersionNumber', False, pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 8, 6), )

    
    VersionNumber = property(__VersionNumber.value, __VersionNumber.set, None, None)

    
    # Element Description uses Python identifier Description
    __Description = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'Description'), 'Description', '__AbsentNamespace0_AresysXmlDocType_Description', False, pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 9, 6), )

    
    Description = property(__Description.value, __Description.set, None, None)

    
    # Element Channel uses Python identifier Channel
    __Channel = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'Channel'), 'Channel', '__AbsentNamespace0_AresysXmlDocType_Channel', True, pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 11, 8), )

    
    Channel = property(__Channel.value, __Channel.set, None, None)

    _ElementMap.update({
        __NumberOfChannels.name() : __NumberOfChannels,
        __VersionNumber.name() : __VersionNumber,
        __Description.name() : __Description,
        __Channel.name() : __Channel
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.AresysXmlDocType = AresysXmlDocType
Namespace.addCategoryObject('typeBinding', 'AresysXmlDocType', AresysXmlDocType)


# Complex type ChannelType with content type EMPTY
class ChannelType (_ImportedBinding_pyxb_metadata_types.TreeElementBaseType):
    """Complex type ChannelType with content type EMPTY"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_EMPTY
    _Abstract = False
    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'ChannelType')
    _XSDLocation = pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 45, 2)
    _ElementMap = _ImportedBinding_pyxb_metadata_types.TreeElementBaseType._ElementMap.copy()
    _AttributeMap = _ImportedBinding_pyxb_metadata_types.TreeElementBaseType._AttributeMap.copy()
    # Base type is _ImportedBinding_pyxb_metadata_types.TreeElementBaseType
    
    # Attribute ContentID uses Python identifier ContentID
    __ContentID = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'ContentID'), 'ContentID', '__AbsentNamespace0_ChannelType_ContentID', pyxb.binding.datatypes.string)
    __ContentID._DeclarationLocation = pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 48, 8)
    __ContentID._UseLocation = pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 48, 8)
    
    ContentID = property(__ContentID.value, __ContentID.set, None, None)

    
    # Attribute Number inherited from {aresysTypes}TreeElementBaseType
    
    # Attribute Total inherited from {aresysTypes}TreeElementBaseType
    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __ContentID.name() : __ContentID
    })
_module_typeBindings.ChannelType = ChannelType
Namespace.addCategoryObject('typeBinding', 'ChannelType', ChannelType)


# Complex type [anonymous] with content type ELEMENT_ONLY
class CTD_ANON (ChannelType):
    """Complex type [anonymous] with content type ELEMENT_ONLY"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 12, 10)
    _ElementMap = ChannelType._ElementMap.copy()
    _AttributeMap = ChannelType._AttributeMap.copy()
    # Base type is ChannelType
    
    # Element RasterInfo uses Python identifier RasterInfo
    __RasterInfo = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'RasterInfo'), 'RasterInfo', '__AbsentNamespace0_CTD_ANON_RasterInfo', False, pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 16, 18), )

    
    RasterInfo = property(__RasterInfo.value, __RasterInfo.set, None, None)

    
    # Element ROI uses Python identifier ROI
    __ROI = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'ROI'), 'ROI', '__AbsentNamespace0_CTD_ANON_ROI', False, pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 17, 18), )

    
    ROI = property(__ROI.value, __ROI.set, None, None)

    
    # Element DataSetInfo uses Python identifier DataSetInfo
    __DataSetInfo = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'DataSetInfo'), 'DataSetInfo', '__AbsentNamespace0_CTD_ANON_DataSetInfo', False, pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 18, 18), )

    
    DataSetInfo = property(__DataSetInfo.value, __DataSetInfo.set, None, None)

    
    # Element SwathInfo uses Python identifier SwathInfo
    __SwathInfo = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'SwathInfo'), 'SwathInfo', '__AbsentNamespace0_CTD_ANON_SwathInfo', False, pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 19, 18), )

    
    SwathInfo = property(__SwathInfo.value, __SwathInfo.set, None, None)

    
    # Element SamplingConstants uses Python identifier SamplingConstants
    __SamplingConstants = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'SamplingConstants'), 'SamplingConstants', '__AbsentNamespace0_CTD_ANON_SamplingConstants', False, pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 20, 18), )

    
    SamplingConstants = property(__SamplingConstants.value, __SamplingConstants.set, None, None)

    
    # Element AcquisitionTimeLine uses Python identifier AcquisitionTimeLine
    __AcquisitionTimeLine = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'AcquisitionTimeLine'), 'AcquisitionTimeLine', '__AbsentNamespace0_CTD_ANON_AcquisitionTimeLine', False, pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 21, 18), )

    
    AcquisitionTimeLine = property(__AcquisitionTimeLine.value, __AcquisitionTimeLine.set, None, None)

    
    # Element DataStatistics uses Python identifier DataStatistics
    __DataStatistics = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'DataStatistics'), 'DataStatistics', '__AbsentNamespace0_CTD_ANON_DataStatistics', False, pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 22, 18), )

    
    DataStatistics = property(__DataStatistics.value, __DataStatistics.set, None, None)

    
    # Element BurstInfo uses Python identifier BurstInfo
    __BurstInfo = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'BurstInfo'), 'BurstInfo', '__AbsentNamespace0_CTD_ANON_BurstInfo', False, pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 23, 18), )

    
    BurstInfo = property(__BurstInfo.value, __BurstInfo.set, None, None)

    
    # Element StateVectorData uses Python identifier StateVectorData
    __StateVectorData = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'StateVectorData'), 'StateVectorData', '__AbsentNamespace0_CTD_ANON_StateVectorData', False, pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 24, 18), )

    
    StateVectorData = property(__StateVectorData.value, __StateVectorData.set, None, None)

    
    # Element DopplerCentroid uses Python identifier DopplerCentroid
    __DopplerCentroid = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'DopplerCentroid'), 'DopplerCentroid', '__AbsentNamespace0_CTD_ANON_DopplerCentroid', True, pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 25, 18), )

    
    DopplerCentroid = property(__DopplerCentroid.value, __DopplerCentroid.set, None, None)

    
    # Element DopplerRate uses Python identifier DopplerRate
    __DopplerRate = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'DopplerRate'), 'DopplerRate', '__AbsentNamespace0_CTD_ANON_DopplerRate', True, pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 26, 18), )

    
    DopplerRate = property(__DopplerRate.value, __DopplerRate.set, None, None)

    
    # Element TopsAzimuthModulationRate uses Python identifier TopsAzimuthModulationRate
    __TopsAzimuthModulationRate = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'TopsAzimuthModulationRate'), 'TopsAzimuthModulationRate', '__AbsentNamespace0_CTD_ANON_TopsAzimuthModulationRate', True, pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 27, 18), )

    
    TopsAzimuthModulationRate = property(__TopsAzimuthModulationRate.value, __TopsAzimuthModulationRate.set, None, None)

    
    # Element SlantToGround uses Python identifier SlantToGround
    __SlantToGround = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'SlantToGround'), 'SlantToGround', '__AbsentNamespace0_CTD_ANON_SlantToGround', True, pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 28, 18), )

    
    SlantToGround = property(__SlantToGround.value, __SlantToGround.set, None, None)

    
    # Element GroundToSlant uses Python identifier GroundToSlant
    __GroundToSlant = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'GroundToSlant'), 'GroundToSlant', '__AbsentNamespace0_CTD_ANON_GroundToSlant', True, pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 29, 18), )

    
    GroundToSlant = property(__GroundToSlant.value, __GroundToSlant.set, None, None)

    
    # Element SlantToIncidence uses Python identifier SlantToIncidence
    __SlantToIncidence = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'SlantToIncidence'), 'SlantToIncidence', '__AbsentNamespace0_CTD_ANON_SlantToIncidence', True, pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 30, 18), )

    
    SlantToIncidence = property(__SlantToIncidence.value, __SlantToIncidence.set, None, None)

    
    # Element SlantToElevation uses Python identifier SlantToElevation
    __SlantToElevation = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'SlantToElevation'), 'SlantToElevation', '__AbsentNamespace0_CTD_ANON_SlantToElevation', True, pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 31, 18), )

    
    SlantToElevation = property(__SlantToElevation.value, __SlantToElevation.set, None, None)

    
    # Element AttitudeInfo uses Python identifier AttitudeInfo
    __AttitudeInfo = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'AttitudeInfo'), 'AttitudeInfo', '__AbsentNamespace0_CTD_ANON_AttitudeInfo', False, pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 32, 18), )

    
    AttitudeInfo = property(__AttitudeInfo.value, __AttitudeInfo.set, None, None)

    
    # Element GroundCornerPoints uses Python identifier GroundCornerPoints
    __GroundCornerPoints = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'GroundCornerPoints'), 'GroundCornerPoints', '__AbsentNamespace0_CTD_ANON_GroundCornerPoints', False, pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 33, 18), )

    
    GroundCornerPoints = property(__GroundCornerPoints.value, __GroundCornerPoints.set, None, None)

    
    # Element Pulse uses Python identifier Pulse
    __Pulse = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'Pulse'), 'Pulse', '__AbsentNamespace0_CTD_ANON_Pulse', False, pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 34, 18), )

    
    Pulse = property(__Pulse.value, __Pulse.set, None, None)

    
    # Element CoregPoly uses Python identifier CoregPoly
    __CoregPoly = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'CoregPoly'), 'CoregPoly', '__AbsentNamespace0_CTD_ANON_CoregPoly', True, pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 35, 18), )

    
    CoregPoly = property(__CoregPoly.value, __CoregPoly.set, None, None)

    
    # Element AntennaInfo uses Python identifier AntennaInfo
    __AntennaInfo = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'AntennaInfo'), 'AntennaInfo', '__AbsentNamespace0_CTD_ANON_AntennaInfo', False, pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 36, 18), )

    
    AntennaInfo = property(__AntennaInfo.value, __AntennaInfo.set, None, None)

    
    # Attribute ContentID inherited from ChannelType
    
    # Attribute Number inherited from {aresysTypes}TreeElementBaseType
    
    # Attribute Total inherited from {aresysTypes}TreeElementBaseType
    _ElementMap.update({
        __RasterInfo.name() : __RasterInfo,
        __ROI.name() : __ROI,
        __DataSetInfo.name() : __DataSetInfo,
        __SwathInfo.name() : __SwathInfo,
        __SamplingConstants.name() : __SamplingConstants,
        __AcquisitionTimeLine.name() : __AcquisitionTimeLine,
        __DataStatistics.name() : __DataStatistics,
        __BurstInfo.name() : __BurstInfo,
        __StateVectorData.name() : __StateVectorData,
        __DopplerCentroid.name() : __DopplerCentroid,
        __DopplerRate.name() : __DopplerRate,
        __TopsAzimuthModulationRate.name() : __TopsAzimuthModulationRate,
        __SlantToGround.name() : __SlantToGround,
        __GroundToSlant.name() : __GroundToSlant,
        __SlantToIncidence.name() : __SlantToIncidence,
        __SlantToElevation.name() : __SlantToElevation,
        __AttitudeInfo.name() : __AttitudeInfo,
        __GroundCornerPoints.name() : __GroundCornerPoints,
        __Pulse.name() : __Pulse,
        __CoregPoly.name() : __CoregPoly,
        __AntennaInfo.name() : __AntennaInfo
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.CTD_ANON = CTD_ANON


AresysXmlDoc = pyxb.binding.basis.element(pyxb.namespace.ExpandedName(Namespace, 'AresysXmlDoc'), AresysXmlDocType, location=pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 4, 2))
Namespace.addCategoryObject('elementBinding', AresysXmlDoc.name().localName(), AresysXmlDoc)



AresysXmlDocType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'NumberOfChannels'), pyxb.binding.datatypes.unsignedInt, scope=AresysXmlDocType, location=pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 7, 6)))

AresysXmlDocType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'VersionNumber'), pyxb.binding.datatypes.double, scope=AresysXmlDocType, location=pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 8, 6)))

AresysXmlDocType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'Description'), pyxb.binding.datatypes.string, scope=AresysXmlDocType, location=pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 9, 6)))

AresysXmlDocType._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'Channel'), CTD_ANON, scope=AresysXmlDocType, location=pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 11, 8)))

def _BuildAutomaton ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton
    del _BuildAutomaton
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=0, max=None, metadata=pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 10, 6))
    counters.add(cc_0)
    states = []
    final_update = None
    symbol = pyxb.binding.content.ElementUse(AresysXmlDocType._UseForTag(pyxb.namespace.ExpandedName(None, 'NumberOfChannels')), pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 7, 6))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    final_update = None
    symbol = pyxb.binding.content.ElementUse(AresysXmlDocType._UseForTag(pyxb.namespace.ExpandedName(None, 'VersionNumber')), pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 8, 6))
    st_1 = fac.State(symbol, is_initial=False, final_update=final_update, is_unordered_catenation=False)
    states.append(st_1)
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(AresysXmlDocType._UseForTag(pyxb.namespace.ExpandedName(None, 'Description')), pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 9, 6))
    st_2 = fac.State(symbol, is_initial=False, final_update=final_update, is_unordered_catenation=False)
    states.append(st_2)
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_0, False))
    symbol = pyxb.binding.content.ElementUse(AresysXmlDocType._UseForTag(pyxb.namespace.ExpandedName(None, 'Channel')), pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 11, 8))
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
    transitions.append(fac.Transition(st_3, [
        fac.UpdateInstruction(cc_0, True) ]))
    st_3._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=None)
AresysXmlDocType._Automaton = _BuildAutomaton()




CTD_ANON._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'RasterInfo'), _ImportedBinding_pyxb_metadata_types.RasterInfoType, scope=CTD_ANON, location=pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 16, 18)))

CTD_ANON._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'ROI'), _ImportedBinding_pyxb_metadata_types.ROIType, scope=CTD_ANON, location=pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 17, 18)))

CTD_ANON._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'DataSetInfo'), _ImportedBinding_pyxb_metadata_types.DataSetInfoType, scope=CTD_ANON, location=pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 18, 18)))

CTD_ANON._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'SwathInfo'), _ImportedBinding_pyxb_metadata_types.SwathInfoType, scope=CTD_ANON, location=pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 19, 18)))

CTD_ANON._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'SamplingConstants'), _ImportedBinding_pyxb_metadata_types.SamplingConstantsType, scope=CTD_ANON, location=pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 20, 18)))

CTD_ANON._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'AcquisitionTimeLine'), _ImportedBinding_pyxb_metadata_types.AcquisitionTimelineType, scope=CTD_ANON, location=pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 21, 18)))

CTD_ANON._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'DataStatistics'), _ImportedBinding_pyxb_metadata_types.DataStatisticsType, scope=CTD_ANON, location=pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 22, 18)))

CTD_ANON._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'BurstInfo'), _ImportedBinding_pyxb_metadata_types.BurstInfoType, scope=CTD_ANON, location=pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 23, 18)))

CTD_ANON._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'StateVectorData'), _ImportedBinding_pyxb_metadata_types.StateVectorDataType, scope=CTD_ANON, location=pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 24, 18)))

CTD_ANON._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'DopplerCentroid'), _ImportedBinding_pyxb_metadata_types.polyType, scope=CTD_ANON, location=pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 25, 18)))

CTD_ANON._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'DopplerRate'), _ImportedBinding_pyxb_metadata_types.polyType, scope=CTD_ANON, location=pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 26, 18)))

CTD_ANON._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'TopsAzimuthModulationRate'), _ImportedBinding_pyxb_metadata_types.polyType, scope=CTD_ANON, location=pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 27, 18)))

CTD_ANON._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'SlantToGround'), _ImportedBinding_pyxb_metadata_types.polyType, scope=CTD_ANON, location=pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 28, 18)))

CTD_ANON._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'GroundToSlant'), _ImportedBinding_pyxb_metadata_types.polyType, scope=CTD_ANON, location=pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 29, 18)))

CTD_ANON._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'SlantToIncidence'), _ImportedBinding_pyxb_metadata_types.polyType, scope=CTD_ANON, location=pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 30, 18)))

CTD_ANON._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'SlantToElevation'), _ImportedBinding_pyxb_metadata_types.polyType, scope=CTD_ANON, location=pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 31, 18)))

CTD_ANON._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'AttitudeInfo'), _ImportedBinding_pyxb_metadata_types.AttitudeInfoType, scope=CTD_ANON, location=pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 32, 18)))

CTD_ANON._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'GroundCornerPoints'), _ImportedBinding_pyxb_metadata_types.GroundCornersPointsType, scope=CTD_ANON, location=pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 33, 18)))

CTD_ANON._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'Pulse'), _ImportedBinding_pyxb_metadata_types.PulseType, scope=CTD_ANON, location=pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 34, 18)))

CTD_ANON._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'CoregPoly'), _ImportedBinding_pyxb_metadata_types.polyCoregType, scope=CTD_ANON, location=pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 35, 18)))

CTD_ANON._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'AntennaInfo'), _ImportedBinding_pyxb_metadata_types.AntennaInfoType, scope=CTD_ANON, location=pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 36, 18)))

def _BuildAutomaton_ ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_
    del _BuildAutomaton_
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 16, 18))
    counters.add(cc_0)
    cc_1 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 17, 18))
    counters.add(cc_1)
    cc_2 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 18, 18))
    counters.add(cc_2)
    cc_3 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 19, 18))
    counters.add(cc_3)
    cc_4 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 20, 18))
    counters.add(cc_4)
    cc_5 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 21, 18))
    counters.add(cc_5)
    cc_6 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 22, 18))
    counters.add(cc_6)
    cc_7 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 23, 18))
    counters.add(cc_7)
    cc_8 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 24, 18))
    counters.add(cc_8)
    cc_9 = fac.CounterCondition(min=0, max=None, metadata=pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 25, 18))
    counters.add(cc_9)
    cc_10 = fac.CounterCondition(min=0, max=None, metadata=pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 26, 18))
    counters.add(cc_10)
    cc_11 = fac.CounterCondition(min=0, max=None, metadata=pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 27, 18))
    counters.add(cc_11)
    cc_12 = fac.CounterCondition(min=0, max=None, metadata=pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 28, 18))
    counters.add(cc_12)
    cc_13 = fac.CounterCondition(min=0, max=None, metadata=pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 29, 18))
    counters.add(cc_13)
    cc_14 = fac.CounterCondition(min=0, max=None, metadata=pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 30, 18))
    counters.add(cc_14)
    cc_15 = fac.CounterCondition(min=0, max=None, metadata=pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 31, 18))
    counters.add(cc_15)
    cc_16 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 32, 18))
    counters.add(cc_16)
    cc_17 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 33, 18))
    counters.add(cc_17)
    cc_18 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 34, 18))
    counters.add(cc_18)
    cc_19 = fac.CounterCondition(min=0, max=None, metadata=pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 35, 18))
    counters.add(cc_19)
    cc_20 = fac.CounterCondition(min=0, max=1, metadata=pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 36, 18))
    counters.add(cc_20)
    states = []
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_0, False))
    symbol = pyxb.binding.content.ElementUse(CTD_ANON._UseForTag(pyxb.namespace.ExpandedName(None, 'RasterInfo')), pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 16, 18))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_1, False))
    symbol = pyxb.binding.content.ElementUse(CTD_ANON._UseForTag(pyxb.namespace.ExpandedName(None, 'ROI')), pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 17, 18))
    st_1 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_1)
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_2, False))
    symbol = pyxb.binding.content.ElementUse(CTD_ANON._UseForTag(pyxb.namespace.ExpandedName(None, 'DataSetInfo')), pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 18, 18))
    st_2 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_2)
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_3, False))
    symbol = pyxb.binding.content.ElementUse(CTD_ANON._UseForTag(pyxb.namespace.ExpandedName(None, 'SwathInfo')), pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 19, 18))
    st_3 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_3)
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_4, False))
    symbol = pyxb.binding.content.ElementUse(CTD_ANON._UseForTag(pyxb.namespace.ExpandedName(None, 'SamplingConstants')), pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 20, 18))
    st_4 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_4)
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_5, False))
    symbol = pyxb.binding.content.ElementUse(CTD_ANON._UseForTag(pyxb.namespace.ExpandedName(None, 'AcquisitionTimeLine')), pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 21, 18))
    st_5 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_5)
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_6, False))
    symbol = pyxb.binding.content.ElementUse(CTD_ANON._UseForTag(pyxb.namespace.ExpandedName(None, 'DataStatistics')), pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 22, 18))
    st_6 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_6)
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_7, False))
    symbol = pyxb.binding.content.ElementUse(CTD_ANON._UseForTag(pyxb.namespace.ExpandedName(None, 'BurstInfo')), pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 23, 18))
    st_7 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_7)
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_8, False))
    symbol = pyxb.binding.content.ElementUse(CTD_ANON._UseForTag(pyxb.namespace.ExpandedName(None, 'StateVectorData')), pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 24, 18))
    st_8 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_8)
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_9, False))
    symbol = pyxb.binding.content.ElementUse(CTD_ANON._UseForTag(pyxb.namespace.ExpandedName(None, 'DopplerCentroid')), pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 25, 18))
    st_9 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_9)
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_10, False))
    symbol = pyxb.binding.content.ElementUse(CTD_ANON._UseForTag(pyxb.namespace.ExpandedName(None, 'DopplerRate')), pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 26, 18))
    st_10 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_10)
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_11, False))
    symbol = pyxb.binding.content.ElementUse(CTD_ANON._UseForTag(pyxb.namespace.ExpandedName(None, 'TopsAzimuthModulationRate')), pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 27, 18))
    st_11 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_11)
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_12, False))
    symbol = pyxb.binding.content.ElementUse(CTD_ANON._UseForTag(pyxb.namespace.ExpandedName(None, 'SlantToGround')), pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 28, 18))
    st_12 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_12)
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_13, False))
    symbol = pyxb.binding.content.ElementUse(CTD_ANON._UseForTag(pyxb.namespace.ExpandedName(None, 'GroundToSlant')), pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 29, 18))
    st_13 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_13)
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_14, False))
    symbol = pyxb.binding.content.ElementUse(CTD_ANON._UseForTag(pyxb.namespace.ExpandedName(None, 'SlantToIncidence')), pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 30, 18))
    st_14 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_14)
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_15, False))
    symbol = pyxb.binding.content.ElementUse(CTD_ANON._UseForTag(pyxb.namespace.ExpandedName(None, 'SlantToElevation')), pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 31, 18))
    st_15 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_15)
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_16, False))
    symbol = pyxb.binding.content.ElementUse(CTD_ANON._UseForTag(pyxb.namespace.ExpandedName(None, 'AttitudeInfo')), pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 32, 18))
    st_16 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_16)
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_17, False))
    symbol = pyxb.binding.content.ElementUse(CTD_ANON._UseForTag(pyxb.namespace.ExpandedName(None, 'GroundCornerPoints')), pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 33, 18))
    st_17 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_17)
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_18, False))
    symbol = pyxb.binding.content.ElementUse(CTD_ANON._UseForTag(pyxb.namespace.ExpandedName(None, 'Pulse')), pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 34, 18))
    st_18 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_18)
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_19, False))
    symbol = pyxb.binding.content.ElementUse(CTD_ANON._UseForTag(pyxb.namespace.ExpandedName(None, 'CoregPoly')), pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 35, 18))
    st_19 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_19)
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_20, False))
    symbol = pyxb.binding.content.ElementUse(CTD_ANON._UseForTag(pyxb.namespace.ExpandedName(None, 'AntennaInfo')), pyxb.utils.utility.Location('aresys_generic_metadata.xsd', 36, 18))
    st_20 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_20)
    transitions = []
    transitions.append(fac.Transition(st_0, [
        fac.UpdateInstruction(cc_0, True) ]))
    transitions.append(fac.Transition(st_1, [
        fac.UpdateInstruction(cc_0, False) ]))
    transitions.append(fac.Transition(st_2, [
        fac.UpdateInstruction(cc_0, False) ]))
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
    transitions.append(fac.Transition(st_8, [
        fac.UpdateInstruction(cc_0, False) ]))
    transitions.append(fac.Transition(st_9, [
        fac.UpdateInstruction(cc_0, False) ]))
    transitions.append(fac.Transition(st_10, [
        fac.UpdateInstruction(cc_0, False) ]))
    transitions.append(fac.Transition(st_11, [
        fac.UpdateInstruction(cc_0, False) ]))
    transitions.append(fac.Transition(st_12, [
        fac.UpdateInstruction(cc_0, False) ]))
    transitions.append(fac.Transition(st_13, [
        fac.UpdateInstruction(cc_0, False) ]))
    transitions.append(fac.Transition(st_14, [
        fac.UpdateInstruction(cc_0, False) ]))
    transitions.append(fac.Transition(st_15, [
        fac.UpdateInstruction(cc_0, False) ]))
    transitions.append(fac.Transition(st_16, [
        fac.UpdateInstruction(cc_0, False) ]))
    transitions.append(fac.Transition(st_17, [
        fac.UpdateInstruction(cc_0, False) ]))
    transitions.append(fac.Transition(st_18, [
        fac.UpdateInstruction(cc_0, False) ]))
    transitions.append(fac.Transition(st_19, [
        fac.UpdateInstruction(cc_0, False) ]))
    transitions.append(fac.Transition(st_20, [
        fac.UpdateInstruction(cc_0, False) ]))
    st_0._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_1, [
        fac.UpdateInstruction(cc_1, True) ]))
    transitions.append(fac.Transition(st_2, [
        fac.UpdateInstruction(cc_1, False) ]))
    transitions.append(fac.Transition(st_3, [
        fac.UpdateInstruction(cc_1, False) ]))
    transitions.append(fac.Transition(st_4, [
        fac.UpdateInstruction(cc_1, False) ]))
    transitions.append(fac.Transition(st_5, [
        fac.UpdateInstruction(cc_1, False) ]))
    transitions.append(fac.Transition(st_6, [
        fac.UpdateInstruction(cc_1, False) ]))
    transitions.append(fac.Transition(st_7, [
        fac.UpdateInstruction(cc_1, False) ]))
    transitions.append(fac.Transition(st_8, [
        fac.UpdateInstruction(cc_1, False) ]))
    transitions.append(fac.Transition(st_9, [
        fac.UpdateInstruction(cc_1, False) ]))
    transitions.append(fac.Transition(st_10, [
        fac.UpdateInstruction(cc_1, False) ]))
    transitions.append(fac.Transition(st_11, [
        fac.UpdateInstruction(cc_1, False) ]))
    transitions.append(fac.Transition(st_12, [
        fac.UpdateInstruction(cc_1, False) ]))
    transitions.append(fac.Transition(st_13, [
        fac.UpdateInstruction(cc_1, False) ]))
    transitions.append(fac.Transition(st_14, [
        fac.UpdateInstruction(cc_1, False) ]))
    transitions.append(fac.Transition(st_15, [
        fac.UpdateInstruction(cc_1, False) ]))
    transitions.append(fac.Transition(st_16, [
        fac.UpdateInstruction(cc_1, False) ]))
    transitions.append(fac.Transition(st_17, [
        fac.UpdateInstruction(cc_1, False) ]))
    transitions.append(fac.Transition(st_18, [
        fac.UpdateInstruction(cc_1, False) ]))
    transitions.append(fac.Transition(st_19, [
        fac.UpdateInstruction(cc_1, False) ]))
    transitions.append(fac.Transition(st_20, [
        fac.UpdateInstruction(cc_1, False) ]))
    st_1._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_2, [
        fac.UpdateInstruction(cc_2, True) ]))
    transitions.append(fac.Transition(st_3, [
        fac.UpdateInstruction(cc_2, False) ]))
    transitions.append(fac.Transition(st_4, [
        fac.UpdateInstruction(cc_2, False) ]))
    transitions.append(fac.Transition(st_5, [
        fac.UpdateInstruction(cc_2, False) ]))
    transitions.append(fac.Transition(st_6, [
        fac.UpdateInstruction(cc_2, False) ]))
    transitions.append(fac.Transition(st_7, [
        fac.UpdateInstruction(cc_2, False) ]))
    transitions.append(fac.Transition(st_8, [
        fac.UpdateInstruction(cc_2, False) ]))
    transitions.append(fac.Transition(st_9, [
        fac.UpdateInstruction(cc_2, False) ]))
    transitions.append(fac.Transition(st_10, [
        fac.UpdateInstruction(cc_2, False) ]))
    transitions.append(fac.Transition(st_11, [
        fac.UpdateInstruction(cc_2, False) ]))
    transitions.append(fac.Transition(st_12, [
        fac.UpdateInstruction(cc_2, False) ]))
    transitions.append(fac.Transition(st_13, [
        fac.UpdateInstruction(cc_2, False) ]))
    transitions.append(fac.Transition(st_14, [
        fac.UpdateInstruction(cc_2, False) ]))
    transitions.append(fac.Transition(st_15, [
        fac.UpdateInstruction(cc_2, False) ]))
    transitions.append(fac.Transition(st_16, [
        fac.UpdateInstruction(cc_2, False) ]))
    transitions.append(fac.Transition(st_17, [
        fac.UpdateInstruction(cc_2, False) ]))
    transitions.append(fac.Transition(st_18, [
        fac.UpdateInstruction(cc_2, False) ]))
    transitions.append(fac.Transition(st_19, [
        fac.UpdateInstruction(cc_2, False) ]))
    transitions.append(fac.Transition(st_20, [
        fac.UpdateInstruction(cc_2, False) ]))
    st_2._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_3, [
        fac.UpdateInstruction(cc_3, True) ]))
    transitions.append(fac.Transition(st_4, [
        fac.UpdateInstruction(cc_3, False) ]))
    transitions.append(fac.Transition(st_5, [
        fac.UpdateInstruction(cc_3, False) ]))
    transitions.append(fac.Transition(st_6, [
        fac.UpdateInstruction(cc_3, False) ]))
    transitions.append(fac.Transition(st_7, [
        fac.UpdateInstruction(cc_3, False) ]))
    transitions.append(fac.Transition(st_8, [
        fac.UpdateInstruction(cc_3, False) ]))
    transitions.append(fac.Transition(st_9, [
        fac.UpdateInstruction(cc_3, False) ]))
    transitions.append(fac.Transition(st_10, [
        fac.UpdateInstruction(cc_3, False) ]))
    transitions.append(fac.Transition(st_11, [
        fac.UpdateInstruction(cc_3, False) ]))
    transitions.append(fac.Transition(st_12, [
        fac.UpdateInstruction(cc_3, False) ]))
    transitions.append(fac.Transition(st_13, [
        fac.UpdateInstruction(cc_3, False) ]))
    transitions.append(fac.Transition(st_14, [
        fac.UpdateInstruction(cc_3, False) ]))
    transitions.append(fac.Transition(st_15, [
        fac.UpdateInstruction(cc_3, False) ]))
    transitions.append(fac.Transition(st_16, [
        fac.UpdateInstruction(cc_3, False) ]))
    transitions.append(fac.Transition(st_17, [
        fac.UpdateInstruction(cc_3, False) ]))
    transitions.append(fac.Transition(st_18, [
        fac.UpdateInstruction(cc_3, False) ]))
    transitions.append(fac.Transition(st_19, [
        fac.UpdateInstruction(cc_3, False) ]))
    transitions.append(fac.Transition(st_20, [
        fac.UpdateInstruction(cc_3, False) ]))
    st_3._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_4, [
        fac.UpdateInstruction(cc_4, True) ]))
    transitions.append(fac.Transition(st_5, [
        fac.UpdateInstruction(cc_4, False) ]))
    transitions.append(fac.Transition(st_6, [
        fac.UpdateInstruction(cc_4, False) ]))
    transitions.append(fac.Transition(st_7, [
        fac.UpdateInstruction(cc_4, False) ]))
    transitions.append(fac.Transition(st_8, [
        fac.UpdateInstruction(cc_4, False) ]))
    transitions.append(fac.Transition(st_9, [
        fac.UpdateInstruction(cc_4, False) ]))
    transitions.append(fac.Transition(st_10, [
        fac.UpdateInstruction(cc_4, False) ]))
    transitions.append(fac.Transition(st_11, [
        fac.UpdateInstruction(cc_4, False) ]))
    transitions.append(fac.Transition(st_12, [
        fac.UpdateInstruction(cc_4, False) ]))
    transitions.append(fac.Transition(st_13, [
        fac.UpdateInstruction(cc_4, False) ]))
    transitions.append(fac.Transition(st_14, [
        fac.UpdateInstruction(cc_4, False) ]))
    transitions.append(fac.Transition(st_15, [
        fac.UpdateInstruction(cc_4, False) ]))
    transitions.append(fac.Transition(st_16, [
        fac.UpdateInstruction(cc_4, False) ]))
    transitions.append(fac.Transition(st_17, [
        fac.UpdateInstruction(cc_4, False) ]))
    transitions.append(fac.Transition(st_18, [
        fac.UpdateInstruction(cc_4, False) ]))
    transitions.append(fac.Transition(st_19, [
        fac.UpdateInstruction(cc_4, False) ]))
    transitions.append(fac.Transition(st_20, [
        fac.UpdateInstruction(cc_4, False) ]))
    st_4._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_5, [
        fac.UpdateInstruction(cc_5, True) ]))
    transitions.append(fac.Transition(st_6, [
        fac.UpdateInstruction(cc_5, False) ]))
    transitions.append(fac.Transition(st_7, [
        fac.UpdateInstruction(cc_5, False) ]))
    transitions.append(fac.Transition(st_8, [
        fac.UpdateInstruction(cc_5, False) ]))
    transitions.append(fac.Transition(st_9, [
        fac.UpdateInstruction(cc_5, False) ]))
    transitions.append(fac.Transition(st_10, [
        fac.UpdateInstruction(cc_5, False) ]))
    transitions.append(fac.Transition(st_11, [
        fac.UpdateInstruction(cc_5, False) ]))
    transitions.append(fac.Transition(st_12, [
        fac.UpdateInstruction(cc_5, False) ]))
    transitions.append(fac.Transition(st_13, [
        fac.UpdateInstruction(cc_5, False) ]))
    transitions.append(fac.Transition(st_14, [
        fac.UpdateInstruction(cc_5, False) ]))
    transitions.append(fac.Transition(st_15, [
        fac.UpdateInstruction(cc_5, False) ]))
    transitions.append(fac.Transition(st_16, [
        fac.UpdateInstruction(cc_5, False) ]))
    transitions.append(fac.Transition(st_17, [
        fac.UpdateInstruction(cc_5, False) ]))
    transitions.append(fac.Transition(st_18, [
        fac.UpdateInstruction(cc_5, False) ]))
    transitions.append(fac.Transition(st_19, [
        fac.UpdateInstruction(cc_5, False) ]))
    transitions.append(fac.Transition(st_20, [
        fac.UpdateInstruction(cc_5, False) ]))
    st_5._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_6, [
        fac.UpdateInstruction(cc_6, True) ]))
    transitions.append(fac.Transition(st_7, [
        fac.UpdateInstruction(cc_6, False) ]))
    transitions.append(fac.Transition(st_8, [
        fac.UpdateInstruction(cc_6, False) ]))
    transitions.append(fac.Transition(st_9, [
        fac.UpdateInstruction(cc_6, False) ]))
    transitions.append(fac.Transition(st_10, [
        fac.UpdateInstruction(cc_6, False) ]))
    transitions.append(fac.Transition(st_11, [
        fac.UpdateInstruction(cc_6, False) ]))
    transitions.append(fac.Transition(st_12, [
        fac.UpdateInstruction(cc_6, False) ]))
    transitions.append(fac.Transition(st_13, [
        fac.UpdateInstruction(cc_6, False) ]))
    transitions.append(fac.Transition(st_14, [
        fac.UpdateInstruction(cc_6, False) ]))
    transitions.append(fac.Transition(st_15, [
        fac.UpdateInstruction(cc_6, False) ]))
    transitions.append(fac.Transition(st_16, [
        fac.UpdateInstruction(cc_6, False) ]))
    transitions.append(fac.Transition(st_17, [
        fac.UpdateInstruction(cc_6, False) ]))
    transitions.append(fac.Transition(st_18, [
        fac.UpdateInstruction(cc_6, False) ]))
    transitions.append(fac.Transition(st_19, [
        fac.UpdateInstruction(cc_6, False) ]))
    transitions.append(fac.Transition(st_20, [
        fac.UpdateInstruction(cc_6, False) ]))
    st_6._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_7, [
        fac.UpdateInstruction(cc_7, True) ]))
    transitions.append(fac.Transition(st_8, [
        fac.UpdateInstruction(cc_7, False) ]))
    transitions.append(fac.Transition(st_9, [
        fac.UpdateInstruction(cc_7, False) ]))
    transitions.append(fac.Transition(st_10, [
        fac.UpdateInstruction(cc_7, False) ]))
    transitions.append(fac.Transition(st_11, [
        fac.UpdateInstruction(cc_7, False) ]))
    transitions.append(fac.Transition(st_12, [
        fac.UpdateInstruction(cc_7, False) ]))
    transitions.append(fac.Transition(st_13, [
        fac.UpdateInstruction(cc_7, False) ]))
    transitions.append(fac.Transition(st_14, [
        fac.UpdateInstruction(cc_7, False) ]))
    transitions.append(fac.Transition(st_15, [
        fac.UpdateInstruction(cc_7, False) ]))
    transitions.append(fac.Transition(st_16, [
        fac.UpdateInstruction(cc_7, False) ]))
    transitions.append(fac.Transition(st_17, [
        fac.UpdateInstruction(cc_7, False) ]))
    transitions.append(fac.Transition(st_18, [
        fac.UpdateInstruction(cc_7, False) ]))
    transitions.append(fac.Transition(st_19, [
        fac.UpdateInstruction(cc_7, False) ]))
    transitions.append(fac.Transition(st_20, [
        fac.UpdateInstruction(cc_7, False) ]))
    st_7._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_8, [
        fac.UpdateInstruction(cc_8, True) ]))
    transitions.append(fac.Transition(st_9, [
        fac.UpdateInstruction(cc_8, False) ]))
    transitions.append(fac.Transition(st_10, [
        fac.UpdateInstruction(cc_8, False) ]))
    transitions.append(fac.Transition(st_11, [
        fac.UpdateInstruction(cc_8, False) ]))
    transitions.append(fac.Transition(st_12, [
        fac.UpdateInstruction(cc_8, False) ]))
    transitions.append(fac.Transition(st_13, [
        fac.UpdateInstruction(cc_8, False) ]))
    transitions.append(fac.Transition(st_14, [
        fac.UpdateInstruction(cc_8, False) ]))
    transitions.append(fac.Transition(st_15, [
        fac.UpdateInstruction(cc_8, False) ]))
    transitions.append(fac.Transition(st_16, [
        fac.UpdateInstruction(cc_8, False) ]))
    transitions.append(fac.Transition(st_17, [
        fac.UpdateInstruction(cc_8, False) ]))
    transitions.append(fac.Transition(st_18, [
        fac.UpdateInstruction(cc_8, False) ]))
    transitions.append(fac.Transition(st_19, [
        fac.UpdateInstruction(cc_8, False) ]))
    transitions.append(fac.Transition(st_20, [
        fac.UpdateInstruction(cc_8, False) ]))
    st_8._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_9, [
        fac.UpdateInstruction(cc_9, True) ]))
    transitions.append(fac.Transition(st_10, [
        fac.UpdateInstruction(cc_9, False) ]))
    transitions.append(fac.Transition(st_11, [
        fac.UpdateInstruction(cc_9, False) ]))
    transitions.append(fac.Transition(st_12, [
        fac.UpdateInstruction(cc_9, False) ]))
    transitions.append(fac.Transition(st_13, [
        fac.UpdateInstruction(cc_9, False) ]))
    transitions.append(fac.Transition(st_14, [
        fac.UpdateInstruction(cc_9, False) ]))
    transitions.append(fac.Transition(st_15, [
        fac.UpdateInstruction(cc_9, False) ]))
    transitions.append(fac.Transition(st_16, [
        fac.UpdateInstruction(cc_9, False) ]))
    transitions.append(fac.Transition(st_17, [
        fac.UpdateInstruction(cc_9, False) ]))
    transitions.append(fac.Transition(st_18, [
        fac.UpdateInstruction(cc_9, False) ]))
    transitions.append(fac.Transition(st_19, [
        fac.UpdateInstruction(cc_9, False) ]))
    transitions.append(fac.Transition(st_20, [
        fac.UpdateInstruction(cc_9, False) ]))
    st_9._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_10, [
        fac.UpdateInstruction(cc_10, True) ]))
    transitions.append(fac.Transition(st_11, [
        fac.UpdateInstruction(cc_10, False) ]))
    transitions.append(fac.Transition(st_12, [
        fac.UpdateInstruction(cc_10, False) ]))
    transitions.append(fac.Transition(st_13, [
        fac.UpdateInstruction(cc_10, False) ]))
    transitions.append(fac.Transition(st_14, [
        fac.UpdateInstruction(cc_10, False) ]))
    transitions.append(fac.Transition(st_15, [
        fac.UpdateInstruction(cc_10, False) ]))
    transitions.append(fac.Transition(st_16, [
        fac.UpdateInstruction(cc_10, False) ]))
    transitions.append(fac.Transition(st_17, [
        fac.UpdateInstruction(cc_10, False) ]))
    transitions.append(fac.Transition(st_18, [
        fac.UpdateInstruction(cc_10, False) ]))
    transitions.append(fac.Transition(st_19, [
        fac.UpdateInstruction(cc_10, False) ]))
    transitions.append(fac.Transition(st_20, [
        fac.UpdateInstruction(cc_10, False) ]))
    st_10._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_11, [
        fac.UpdateInstruction(cc_11, True) ]))
    transitions.append(fac.Transition(st_12, [
        fac.UpdateInstruction(cc_11, False) ]))
    transitions.append(fac.Transition(st_13, [
        fac.UpdateInstruction(cc_11, False) ]))
    transitions.append(fac.Transition(st_14, [
        fac.UpdateInstruction(cc_11, False) ]))
    transitions.append(fac.Transition(st_15, [
        fac.UpdateInstruction(cc_11, False) ]))
    transitions.append(fac.Transition(st_16, [
        fac.UpdateInstruction(cc_11, False) ]))
    transitions.append(fac.Transition(st_17, [
        fac.UpdateInstruction(cc_11, False) ]))
    transitions.append(fac.Transition(st_18, [
        fac.UpdateInstruction(cc_11, False) ]))
    transitions.append(fac.Transition(st_19, [
        fac.UpdateInstruction(cc_11, False) ]))
    transitions.append(fac.Transition(st_20, [
        fac.UpdateInstruction(cc_11, False) ]))
    st_11._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_12, [
        fac.UpdateInstruction(cc_12, True) ]))
    transitions.append(fac.Transition(st_13, [
        fac.UpdateInstruction(cc_12, False) ]))
    transitions.append(fac.Transition(st_14, [
        fac.UpdateInstruction(cc_12, False) ]))
    transitions.append(fac.Transition(st_15, [
        fac.UpdateInstruction(cc_12, False) ]))
    transitions.append(fac.Transition(st_16, [
        fac.UpdateInstruction(cc_12, False) ]))
    transitions.append(fac.Transition(st_17, [
        fac.UpdateInstruction(cc_12, False) ]))
    transitions.append(fac.Transition(st_18, [
        fac.UpdateInstruction(cc_12, False) ]))
    transitions.append(fac.Transition(st_19, [
        fac.UpdateInstruction(cc_12, False) ]))
    transitions.append(fac.Transition(st_20, [
        fac.UpdateInstruction(cc_12, False) ]))
    st_12._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_13, [
        fac.UpdateInstruction(cc_13, True) ]))
    transitions.append(fac.Transition(st_14, [
        fac.UpdateInstruction(cc_13, False) ]))
    transitions.append(fac.Transition(st_15, [
        fac.UpdateInstruction(cc_13, False) ]))
    transitions.append(fac.Transition(st_16, [
        fac.UpdateInstruction(cc_13, False) ]))
    transitions.append(fac.Transition(st_17, [
        fac.UpdateInstruction(cc_13, False) ]))
    transitions.append(fac.Transition(st_18, [
        fac.UpdateInstruction(cc_13, False) ]))
    transitions.append(fac.Transition(st_19, [
        fac.UpdateInstruction(cc_13, False) ]))
    transitions.append(fac.Transition(st_20, [
        fac.UpdateInstruction(cc_13, False) ]))
    st_13._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_14, [
        fac.UpdateInstruction(cc_14, True) ]))
    transitions.append(fac.Transition(st_15, [
        fac.UpdateInstruction(cc_14, False) ]))
    transitions.append(fac.Transition(st_16, [
        fac.UpdateInstruction(cc_14, False) ]))
    transitions.append(fac.Transition(st_17, [
        fac.UpdateInstruction(cc_14, False) ]))
    transitions.append(fac.Transition(st_18, [
        fac.UpdateInstruction(cc_14, False) ]))
    transitions.append(fac.Transition(st_19, [
        fac.UpdateInstruction(cc_14, False) ]))
    transitions.append(fac.Transition(st_20, [
        fac.UpdateInstruction(cc_14, False) ]))
    st_14._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_15, [
        fac.UpdateInstruction(cc_15, True) ]))
    transitions.append(fac.Transition(st_16, [
        fac.UpdateInstruction(cc_15, False) ]))
    transitions.append(fac.Transition(st_17, [
        fac.UpdateInstruction(cc_15, False) ]))
    transitions.append(fac.Transition(st_18, [
        fac.UpdateInstruction(cc_15, False) ]))
    transitions.append(fac.Transition(st_19, [
        fac.UpdateInstruction(cc_15, False) ]))
    transitions.append(fac.Transition(st_20, [
        fac.UpdateInstruction(cc_15, False) ]))
    st_15._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_16, [
        fac.UpdateInstruction(cc_16, True) ]))
    transitions.append(fac.Transition(st_17, [
        fac.UpdateInstruction(cc_16, False) ]))
    transitions.append(fac.Transition(st_18, [
        fac.UpdateInstruction(cc_16, False) ]))
    transitions.append(fac.Transition(st_19, [
        fac.UpdateInstruction(cc_16, False) ]))
    transitions.append(fac.Transition(st_20, [
        fac.UpdateInstruction(cc_16, False) ]))
    st_16._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_17, [
        fac.UpdateInstruction(cc_17, True) ]))
    transitions.append(fac.Transition(st_18, [
        fac.UpdateInstruction(cc_17, False) ]))
    transitions.append(fac.Transition(st_19, [
        fac.UpdateInstruction(cc_17, False) ]))
    transitions.append(fac.Transition(st_20, [
        fac.UpdateInstruction(cc_17, False) ]))
    st_17._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_18, [
        fac.UpdateInstruction(cc_18, True) ]))
    transitions.append(fac.Transition(st_19, [
        fac.UpdateInstruction(cc_18, False) ]))
    transitions.append(fac.Transition(st_20, [
        fac.UpdateInstruction(cc_18, False) ]))
    st_18._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_19, [
        fac.UpdateInstruction(cc_19, True) ]))
    transitions.append(fac.Transition(st_20, [
        fac.UpdateInstruction(cc_19, False) ]))
    st_19._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_20, [
        fac.UpdateInstruction(cc_20, True) ]))
    st_20._set_transitionSet(transitions)
    return fac.Automaton(states, counters, True, containing_state=None)
CTD_ANON._Automaton = _BuildAutomaton_()

