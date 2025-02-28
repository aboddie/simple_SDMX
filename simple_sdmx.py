"""This module provides useful objects for working with SDMX files.

In particular this module provides a base SDMX object useful for reading the
header or getting the SDMX version. In addition specialized objects for DSDs,
codelists, SDMX time series files, and timeseiries are provided.
"""
from __future__ import annotations

__version__ = '0.1'
__author__ = 'Allen Boddie'

import io
import ssl
#import urllib
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

from typing import Dict
from typing import List
from typing import NamedTuple
from typing import Optional
from typing import Tuple


def append_client_id(url: str, odfa_client_id: str) -> str:
    '''Return url with client id appended if the site id from Open Data for
    Africa otherwise returns url without client id.
    '''
    if (url.find('opendataforafrica.org') > 0):
        if url.find('&client_id') == -1:
            return f'{url}&client_id={odfa_client_id}'
    return url


class Structure_Signature(NamedTuple):
    '''This class is used to represent SDMX structures like DSDs, code lists, 
    dataflows, and provision agreements.
    '''
    stype: str
    agencyID: str
    ID: str
    version: str
    
    def generate_url(self) -> str:
        '''Generates url for definition of structue.'''
        ENTRY_URL = {
            "IMF": "https://sdmxcentral.imf.org/ws/public/sdmxapi/rest/",
            "SDMX": "https://registry.sdmx.org/ws/public/sdmxapi/rest/",
            "UNSD": "http://data.un.org/ws/rest/"
            }
    
        ws_endpoint = ENTRY_URL.get(self.agencyID, ENTRY_URL['SDMX'])
        
        if self.version == None:
            url = (f'{ws_endpoint}{self.stype}/{self.agencyID}/{self.ID}/'
                   f'?format=sdmx-2.1&detail=full&references=none')
        else:
            url = (f'{ws_endpoint}{self.stype}/{self.agencyID}/{self.ID}/'
                   f'{self.version}/?format=sdmx-2.1&detail=full&references=none')
        return url


class SDMX():
    '''This is the base class most of this module is based on. The
    SDMX object accepts either a file path or a url and gives access
    to basic SDMX information: version, headers, namespaces.
    '''

    __slots__ = ['uri', 'namespaces', 'header', 'version']

    def __init__(self, uri: str, timeout: Optional[int] = 10) -> None:
        self.uri = uri
        if urllib.parse.urlparse(uri).scheme in ('http', 'https'):
            headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0;'
                                     'Win64; x64) AppleWebKit/537.36'
                                     '(KHTML, like Gecko) Chrome/'
                                     '75.0.3770.142 Safari/537.36',
                       'Accept': 'text/html,application/xhtml+xml,'
                                 'application/xml;q=0.9,image/'
                                 'webp,image/apng,*/*;q=0.8,application'
                                 '/signed-exchange;v=b3'}
            request = urllib.request.Request(uri, headers=headers)
            try:
                xml = urllib.request.urlopen(request, timeout=timeout, context=ssl._create_unverified_context()).read()
            except urllib.error.HTTPError as error:
                raise Exception(f'{type(self).__name__}: {error.code} Error connecting to url provided: {uri}')
            except urllib.error.URLError:
                raise Exception(f'{type(self).__name__}: Error connecting to url provided: {uri}')
        else:
            try:
                with open(uri, 'rb') as content_file:
                    xml = content_file.read()
            except FileNotFoundError:
                raise Exception(f'{type(self).__name__}: Cannot open file provided: {uri}')
        try:
            tree = ET.ElementTree(ET.fromstring(xml))
            root = tree.getroot()
        except ET.ParseError:
            raise Exception(f'{type(self).__name__}: Cannot parse SDMX file: {uri}')
        self.namespaces = dict(
            [node for _, node in ET.iterparse(io.BytesIO(xml),
                                              events=['start-ns'])])
        self.version = self._get_version(self.namespaces)
        self.header = dict()
        try:
            for element in root.find(self._get_ns('Header'), self.namespaces):
                self.header[element.tag.rpartition('}')[2]] = element.text
        except TypeError: # TODO: SDMX v2.0 empty namespace REMOVE clean up
            try:
                for element in root.find(f'{{{self.namespaces[""]}}}Header'):
                    self.header[element.tag.rpartition('}')[2]] = element.text
            except (TypeError, KeyError):
                try:
                    for element in root.find(f'{{{self.namespaces["message"]}}}Header'):
                        self.header[element.tag.rpartition('}')[2]] = element.text
                except (TypeError, KeyError):
                    #Can't get header with this logic
                    pass   
        self._extra_steps(root, timeout)

    def _extra_steps(self, root: ET.Element, timeout: int) -> None:
        '''Sets up subclasses. use this instead of super().__init__ because
        do not want to carry around root which could be big. There is probably 
        a better way to do this.'''

    def __repr__(self) -> str:
        return f'{type(self).__name__}(\'{self.uri}\')'

    def __str__(self) -> str:
        try:
            return f'{type(self).__name__} File: {self.name}.'
        except AttributeError:
            return f'{type(self).__name__} File.'

    def __hash__(self) -> int:
        return hash(self.uri)

    def __eq__(self, other: SDMX) -> bool:
        if isinstance(other, self.__class__):
            return self.uri == other.uri
        return False

    def _get_ns(self, element_name: str) -> str:
        '''prefix namespace if needed, before adding this we would strip
        namespaces from root'''
        ns_dict = {
            'Header': 'http://www.sdmx.org/resources/'
                      'sdmxml/schemas/v2_1/message',
            'DataStructure': 'http://www.sdmx.org/resources/'
                             'sdmxml/schemas/v2_1/structure',
            'DimensionList': 'http://www.sdmx.org/resources/'
                             'sdmxml/schemas/v2_1/structure',
            'AttributeList': 'http://www.sdmx.org/resources/'
                             'sdmxml/schemas/v2_1/structure',
            'Enumeration': 'http://www.sdmx.org/resources/'
                           'sdmxml/schemas/v2_1/structure',
            'Codelist': 'http://www.sdmx.org/resources/'
                        'sdmxml/schemas/v2_1/structure',
            'Code': 'http://www.sdmx.org/resources/sdmxml'
                    '/schemas/v2_1/structure',
            'Name': 'http://www.sdmx.org/resources/sdmxml/'
                    'schemas/v2_1/common',
            'DataStructureComponents': 'http://www.sdmx.org/resources/'
                                       'sdmxml/schemas/v2_1/structure',
            'ProvisionAgreement': 'http://www.sdmx.org/resources/sdmxml/schemas/v2_1/structure',
            'StructureUsage': 'http://www.sdmx.org/resources/sdmxml/schemas/v2_1/structure',
            'DataProvider': 'http://www.sdmx.org/resources/sdmxml/schemas/v2_1/structure',
            'Dataflow': 'http://www.sdmx.org/resources/sdmxml/schemas/v2_1/structure',
            'Structure':'http://www.sdmx.org/resources/sdmxml/schemas/v2_1/structure'
            }
        for ns_prefix, ns_uri in self.namespaces.items():
            if ns_uri == ns_dict.get(element_name,''):
                return f'{ns_prefix}:{element_name}'
        return f'{element_name}'

    @staticmethod
    def _get_version(namespaces: Dict[str, str]) -> str:
        '''Returns SDMX version by looking within the namedspaces.'''
        phrase = 'http://www.sdmx.org/resources/sdmxml/schemas/'
        for uri in namespaces.values():
            if uri.lower().startswith(phrase):
                version = uri[len(phrase):].split('/')[0]
                if version == 'v2_1':
                    return '2.1'
                elif version == 'v2_0':
                    return '2.0'
        return 'UNKNOWN'


class Dataflow(SDMX):
    '''The Dataflow class is based on the SDMX class. It adds a name of the 
    Dataflow as well as a Structure_Signature for the referenced structure.
    Currently supports Dataflow with one and only one referenced structure.
    '''
    
    __slots__ = ['name','structure']
    
    def _extra_steps(self, root: ET.Element, timeout: int) -> None:
        self.name = root.find(f'.//{self._get_ns("Dataflow")}'
                              f'/{self._get_ns("Name")}', self.namespaces
                              ).text
        structure = root.findall(f'.//{self._get_ns("Dataflow")}'
                              f'/{self._get_ns("Structure")}', self.namespaces
                              )
        if len(structure) == 1:
            structure_attrib = structure[0][0].attrib
            agency_id = structure_attrib['agencyID']
            structure_id = structure_attrib['id']
            structure_version = structure_attrib['version']
            structure_type = structure_attrib['class'].lower()
            self.structure = Structure_Signature(stype = structure_type,
                           agencyID = agency_id,
                           ID = structure_id,
                           version = structure_version)
        else:
            raise Exception(f'{self.name}: Invalid Dataflow, reference to multiple strucutres.')
            
class dataproviders(SDMX):
    __slots__ = ['name', 'providers']
    def _extra_steps(self, root: ET.Element, timeout: int) -> None:
        pass #TODO build out data providers information, similar to codelist
             #and use in Provision Agreements to get humand readable.

class ProvisionAgreements(SDMX):
    '''The ProvisionAgreements class is based on the SDMX class. It adds a name and provider 
    of the ProvisionAgreements as well as a Structure_Signature for the referenced structure.
    '''
    
    __slots__ = ['name','structure', 'provider']
    
    def _extra_steps(self, root: ET.Element, timeout: int) -> None:
        self.name = root.find(f'.//{self._get_ns("ProvisionAgreement")}'
                              f'/{self._get_ns("Name")}', self.namespaces
                              ).text
        structure = root.findall(f'.//{self._get_ns("ProvisionAgreement")}'
                                 f'/{self._get_ns("StructureUsage")}', self.namespaces
                                 )
        if len(structure) == 1:
            structure_attrib = structure[0][0].attrib
            agency_id = structure_attrib['agencyID']
            structure_id = structure_attrib['id']
            structure_version = structure_attrib['version']
            structure_type = structure_attrib['class'].lower()
            self.structure = Structure_Signature(stype=structure_type,
                                                 agencyID=agency_id,
                                                 ID=structure_id,
                                                 version=structure_version)
        else:
            raise Exception(f'{self.name}: Invalid ProvisionAgreements, reference to multiple strucutres.')
        provider = root.findall(f'.//{self._get_ns("ProvisionAgreement")}'
                      f'/{self._get_ns("DataProvider")}', self.namespaces
                      )
        if len(provider) == 1:
            self.provider = provider[0].attrib['id'] #TODO get human readable provider https://registry.sdmx.org/ws/public/sdmxapi/rest/dataproviderscheme/SDMX
        else:
            raise Exception(f'{self.name}: Invalid ProvisionAgreements, reference to multiple providers.')             


class _ValidateSeriesWithDSD():
    '''Private class used to mix these functions into DSD and series.'''

    @staticmethod
    def _validate_series(dsd: DSD, series: SDMXTimeseries) -> bool:
        '''Returns true if all dimensions for DSD are present and populated
        with valid values.
        '''
        for dim, code_list in dsd.dimensions.items():
            if dim in series.metadata.keys():
                if series.metadata[dim] not in code_list.codes():
                    return False
            else:
                return False
        return True

    @staticmethod
    def _validate_series_details(dsd: DSD, series: SDMXTimeseries) -> Dict[str, bool]:
        '''Returns a dict with true or false for each dimension on the DSD.'''
        conformingdict = dict()
        for dim, code_list in dsd.dimensions.items():
            if dim in series.metadata.keys():
                conformingdict[dim] = series.metadata[dim] in code_list.codes()
            else:
                conformingdict[dim] = False
        return conformingdict

    @staticmethod
    def _validate_series_dimnension(dsd: DSD, series: SDMXTimeseries, dimension: str) -> bool:
        '''Returns true if series has a valid entery in the
        supplied dimension.
        '''
        try:
            return (series.metadata[dimension] in
                    dsd.dimensions[dimension].codes())
        except KeyError:
            raise Exception(f'Invalid dimension: {dimension}')

    @staticmethod
    def _name_from_metadata(dsd: DSD, series: SDMXTimeseries) -> Dict[str, str]:
        '''Returns a dictonary with names for dimensions and values.'''
        human_readable = dict()
        for dimension_id, code_list in dsd.dimensions.items():
            human_readable[dimension_id] = code_list.name_from_code(
                series.metadata.get(dimension_id, 'Invalid code'))

            # Not possible to resolve dimension id for example
            # counterpart area (the dimension name) in ECOFIN 
            #TODO use concept schme maybe to get this
        return human_readable
          

class DSD(_ValidateSeriesWithDSD, SDMX):
    '''The DSD class is based on the SDMX class. It adds a dict
    of dimension ids and code lists. This class has methods which can be used
    validate a series or return the names for both dimensions and values for
    a given series.
    '''

    __slots__ = ['name', 'dimensions', 'attributes']

    def _extra_steps(self, root: ET.Element, timeout: int) -> None:
        self.name = root.find(f'.//{self._get_ns("DataStructure")}'
                              f'/{self._get_ns("Name")}', self.namespaces
                              ).text
        self.dimensions = dict()
        for dimension in root.find(f'.//{self._get_ns("DataStructure")}/'
                                   f'{self._get_ns("DataStructureComponents")}'
                                   f'/{self._get_ns("DimensionList")}',
                                   self.namespaces):
            for dimension_ref in dimension.findall(
                    f'.//{self._get_ns("Enumeration")}'
                    f'/Ref', self.namespaces):
                cl = self._generate_codelist_signature(dimension_ref)
                # TODO: This isn't going to work for imbeded codelists
                # such as those produced by EcOS. We still have a list
                # Name and indicator list are all that are important.
                # Maybe make a structure to hold them and use CodeList class
                # to populate if URL otherwise populate from DSD XML.
                self.dimensions[dimension.attrib['id']] = (
                    CodeList(cl.generate_url(), timeout=timeout))
        self.attributes = dict()
        for attribute in root.find(f'.//{self._get_ns("DataStructure")}/'
                           f'{self._get_ns("DataStructureComponents")}'
                           f'/{self._get_ns("AttributeList")}',
                           self.namespaces):
            for attribute_ref in attribute.findall(
                    f'.//{self._get_ns("Enumeration")}'
                    f'/Ref', self.namespaces):
                cl = self._generate_codelist_signature(attribute_ref)
                # TODO: This isn't going to work for imbeded codelists
                self.attributes[attribute.attrib['id']] = (
                    CodeList(cl.generate_url(), timeout=timeout))
                # TODO: Capture Attribute type (Series, Dataset, Observation)
        # TODO: Very Low priority, measures

    def __len__(self) -> int:
        return len(self.dimensions)

    @staticmethod
    def _generate_codelist_signature(cl_ref: ET.Element) -> Structure_Signature:
        cl_agency = cl_ref.attrib['agencyID']
        cl_version = cl_ref.attrib['version']
        cl_id = cl_ref.attrib['id']
        return Structure_Signature(stype = 'codelist',
               agencyID = cl_agency,
               ID = cl_id,
               version = cl_version)

    def validate_series(self, series: SDMXTimeseries) -> bool:
        '''Returns true if all dimensions for DSD are present and populated
        with valid values.
        '''
        return self._validate_series(self, series)

    def validate_series_details(self, series: SDMXTimeseries) -> Dict[str, bool]:
        '''Returns a dict with true or false for each dimension.'''
        return self._validate_series_details(self, series)

    def validate_series_dimnension(self, series: SDMXTimeseries, dimension: str) -> bool:
        '''Returns true if series has a valid entery in the given dimension.'''
        return self._validate_series_dimnension(self, series, dimension)

    def name_from_metadata(self, series: SDMXTimeseries) -> Dict[str, str]:
        '''Returns a dictonary with names for dimensions and values.'''
        return self._name_from_metadata(self, series)


class CodeList(SDMX):
    '''The CodeList class is based on the SDMX class. It adds a dict
    of code ids and names which can be used validate a code or return
    the name for a code.
    '''

    __slots__ = ['name', 'indicators']

    def _extra_steps(self, root: ET.Element, timeout: int) -> None:
        self.name = root.find(f'.//{self._get_ns("Codelist")}'
                              f'/{self._get_ns("Name")}',
                              self.namespaces).text
        self.indicators = dict()
        for code in root.findall(f'.//{self._get_ns("Codelist")}'
                                 f'/{self._get_ns("Code")}',
                                 self.namespaces):
            self.indicators[code.attrib['id']] = code.find(
                f'.//{self._get_ns("Name")}', self.namespaces).text

    def __len__(self) -> int:
        return len(self.indicators)

    def codes(self) -> Tuple[str]:
        '''Returns all codes for indicators list.'''
        return tuple([*self.indicators])

    def validate_code(self, code: str) -> bool:
        '''Returns true if the supplied code is in the code list.'''
        return code in self.indicators.keys()

    def name_from_code(self, code: str) -> str:
        '''Returns the name from the code'''
        try:
            return self.indicators[code]
        except KeyError:
            return 'Invalid code'
            #raise Exception(f'Invalid code: {code}')


class SDMXTimeseries(_ValidateSeriesWithDSD):
    '''SDMXTimeseries class allows for basic manipulation for timeseries.'''

    __slots__ = ['metadata', 'frequency', 'scale', 'data', 'invalid_data',
                 'time_period_coverage']

    def __init__(self, xmlseries: ET.Element, version: str, _dsdns: str = '') -> None:
        self.metadata = xmlseries.attrib.copy()
        self.frequency = self.metadata.get('FREQ', 'UNKNOWN')[0]
        self.scale = self.metadata.get('UNIT_MULT', 'UNKNOWN')[0]
        self.data = dict()
        self.invalid_data = dict()
        if version == '2.1':
            searchterm = 'Obs'
        else:
            searchterm = (f'{{{_dsdns}}}Obs')
        for obs in xmlseries.iter(searchterm):
            try:
                self.data[obs.attrib['TIME_PERIOD']] = (
                    float(obs.get('OBS_VALUE')))
            except (ValueError, TypeError):
                self.invalid_data[obs.attrib['TIME_PERIOD']] = (
                    obs.get('OBS_VALUE',''))
        self.time_period_coverage = [*self.data]
        self.time_period_coverage.sort()
        # TODO: Works for A Q M, maybe not other frequencies. If need custom
        # sort fix time_period_coverage (SDMX data file) too.

    def __len__(self) -> int:
        return len(self.data)

    def __repr__(self) -> str:
        return f'SDMX time series: {self.metadata}'

    @property
    def last_observation(self) -> str:
        '''Returns last date with has data.'''
        try:
            return self.time_period_coverage[-1]
        except IndexError:
            return 'N.A.'

    @property
    def has_nonzero_data(self) -> bool:
        '''Returns true if one datapoint contains non-zero data'''
        if len(self.data) == 0:
            return False
        elif min(self.data) == max(self.data) == 0:
            return False
        return True

    @property
    def has_invalid_data(self) -> bool:
        '''Returns if invalid data exists.'''
        return len(self.invalid_data) > 0

    def validate_series(self, dsd: DSD) -> bool:
        '''Returns true if all dimensions for DSD are present and populated
        with valid values.
        '''
        return self._validate_series(dsd, self)

    def validate_series_details(self, dsd: DSD) -> Dict[str, bool]:
        '''Returns a dict with true or false for each dimension.'''
        return self._validate_series_details(dsd, self)

    def validate_series_dimnension(self, dsd: DSD, dimension: str) -> bool:
        '''Returns true if the supplied dimension contains a valid value.'''
        return self._validate_series_dimnension(dsd, self, dimension)

    def name_from_metadata(self, dsd: DSD) -> Dict[str, str]:
        '''Returns a dictonary with names for dimensions and values.'''
        return self._name_from_metadata(dsd, self)

    def metadata_fields(self) -> Tuple[str]:
        '''Returns metadata fields attached to series.'''
        return tuple(self.metadata.keys())

    # TODO: add back +-*/ this will require a new way to init if not XML.


class SDMXDataFile(SDMX):
    '''The SDMXDataFile class is based on the SDMX class. It adds a dsd and
    a list of series (SDMXTimeseries objects) it also includes methods to roll
    up information from the series.
    '''

    __slots__ = ['dsd', 'series']

    def _extra_steps(self, root: ET.Element, timeout: int) -> None:
        uri, structure = self._get_structure_signiture(self.namespaces)
        if structure.stype == 'dataflow':
            url = structure.generate_url()
            self.dsd = Dataflow(url, timeout).structure 
            # TODO: what if the refernce is not a DSD need to check
        elif structure.stype == 'dataprovision':
            url = structure.generate_url()
            #TODO: add support for provision agreements
            raise Exception(f'Provision agreements not support, and generally not publicly accessible: {url}')
        elif structure.stype == 'datastructure':
            self.dsd = structure
        else:
            raise Exception('Could not find DSD declariation in SDMX file.')
        self.series = []
        if self.version == '2.1':
            for series in root.iter('Series'):
                self.series.append(SDMXTimeseries(series,self.version))
        elif self.version == '2.0':
            for series in root.iter(f'{{{uri}}}Series'):
                self.series.append(SDMXTimeseries(series,self.version, uri))

    def __len__(self) -> int:
        return len(self.series)

    @staticmethod
    def _get_structure_signiture(namespaces: Dict[str, str]) -> Structure_Signature:
        phrases = [('dataflow','urn:sdmx:org.sdmx.infomodel.datastructure.dataflow='),
                   ('dataprovision','urn:sdmx:org.sdmx.infomodel.registry.provisionagreement='),
                   ('datastructure','urn:sdmx:org.sdmx.infomodel.datastructure.datastructure='),
                   ('datastructure','urn:sdmx:org.sdmx.infomodel.keyfamily.keyfamily=')]

        for uri in namespaces.values():
            structuretype, _, structuredef = uri.partition('=')
            for phrasetype, phrase in phrases:
                if uri.lower().startswith(phrase):
                    _, _, structuredef = uri.partition('=')
                    agency, _, remaider = structuredef.partition(":")
                    name, _, _ = remaider.partition(':')
                    structureid, _, structure_version = name.partition('(')
                    if len(structure_version) > 0:
                        structure_version = structure_version[:structure_version.rfind(')')]
                        return uri, Structure_Signature(stype = phrasetype,
                               agencyID = agency,
                               ID = structureid,
                               version = structure_version)
                    else:
                        return uri, Structure_Signature(stype = phrasetype,
                               agencyID = agency,
                               ID = structureid,
                               version = None)
        return None, Structure_Signature(stype = False,
                   agencyID = None,
                   ID = None,
                   version = None)

    def generate_dsd(self) -> DSD:
        '''Returns an instance of the DSD used for this SDMX file.'''
        return DSD(self.dsd.generate_url())
        
    def dimensions(self) -> Tuple[str]:
        '''Returns all dimensions used in root.  All series are looped over
        just in case their is an inconsistency in the file.
        '''
        dimension_list = []
        for series in self.series:
            for dimension in series.metadata_fields():
                if dimension not in dimension_list:
                    dimension_list.append(dimension)
        return tuple(dimension_list)

    def unique_dimension_values(self, dimension: str) -> Tuple[str]:
        '''Returns unique values used within the given dimension.'''
        unique_values = set()
        for series in self.series:
            unique_values.add(series.metadata[dimension])
        return tuple(unique_values)

    def _series_counter(self, series_attribute: str,
                        key: Optional[str] = '') -> Dict[str, int]:
        '''Rolls up series attribute into a dict counting occurances.'''
        counting_dict = dict()
        for series in self.series:
            if key != "":
                series_value = getattr(series, series_attribute)[key]
            else:
                series_value = getattr(series, series_attribute)
            if series_value not in counting_dict.keys():
                counting_dict[series_value] = 1
            else:
                counting_dict[series_value] += 1
        return counting_dict

    def frequencies(self) -> Dict[str, int]:
        '''Returns a dictonary of the number of series of each frequency.'''
        return self._series_counter('frequency')

    def last_data_points(self) -> Dict[str, int]:
        '''Returns a dictonary of the number of series having a given date
        as last observation.
        '''
        return self._series_counter('last_observation')

    def scales(self) -> Dict[str, int]:
        '''Returns a dictonary of the number of series of each scale.'''
        return self._series_counter('scale')

    def metadata(self, dimension: str) -> Dict[str, int]:
        '''Returns a dictonary of the number of series of a given value
        in a dimension.
        '''
        try:
            return self._series_counter('metadata', key=dimension)
        except KeyError:
            raise Exception(f'Invalid dimension: {dimension}')

    def time_period_coverage(self, frequency: str) -> List[Tuple[str, int]]:
        '''Returns a dictonary of the numbr of occurances of data for a
        given date.
        '''
        date_list = dict()
        for series in self.series:
            if series.frequency == frequency:
                for date in series.time_period_coverage:
                    if date not in date_list.keys():
                        date_list[date] = 1
                    else:
                        date_list[date] += 1
        return sorted(date_list.items())
