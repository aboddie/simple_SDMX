"""This module provides useful objects for working with SDMX files.

In particular this module provides a base SDMX object useful for reading the
header or getting the SDMX version. In addition specialized objects for DSDs,
codelists, SDMX time series files, and timeseiries are provided.
"""
# from __future__ import annotations
# TODO: with Python 3.7 un comment out and add back typehints for DSD, SDMX,
# and SDMXtimeseries objects _ functions and __eq__.

__version__ = '0.1'
__author__ = 'Allen Boddie'

import io
import urllib
import xml.etree.ElementTree as ET
from typing import Dict
from typing import List
from typing import Optional
from typing import Tuple


# No way to resolve websevice endpoints from SDMX files.
ENTRY_URL = {
        "IMF": "https://sdmxcentral.imf.org/ws/public/sdmxapi/rest/",
        "SDMX": "https://registry.sdmx.org/ws/public/sdmxapi/rest/",
        "ESTAT": "http://ec.europa.eu/eurostat/SDMX/diss-web/rest/",
        "OECD": "http://stats.oecd.org/restsdmx/sdmx.ashx/",
        "UNSD": "http://data.un.org/ws/rest/"
        }


def append_client_id(url: str, client_id: str) -> str:
    """Return url with client id appended if the site was created by Knoema. If
    site was not created by Knoema returns url.
    """
    if (url.find('opendataforafrica.org') > 0) or (url.find('knoema.com') > 0):
        if url.find('&client_id') == -1:
            return f'{url}&client_id={client_id}'
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
                xml = urllib.request.urlopen(request, timeout=timeout).read()
                # Add context=ssl._create_unverified_context()
                # to urlopen if needed
            except urllib.error.HTTPError as error:
                raise Exception(f'{error.code} error connecting'
                                f'to url provided: {uri}')
            except urllib.error.URLError:
                raise Exception(f'Can not connect to url provided: {uri}')
        else:
            try:
                with open(uri, 'rb') as content_file:
                    xml = content_file.read()
            except FileNotFoundError:
                raise Exception(f'Can not open file provided: {uri}')
        try:
            tree = ET.ElementTree(ET.fromstring(xml))
            root = tree.getroot()
        except ET.ParseError:
            raise Exception(f'Can not read SDMX file: {uri}')
        self.namespaces = dict(
            [node for _, node in ET.iterparse(io.BytesIO(xml),
                                              events=['start-ns'])])
        self.version = self._get_version(self.namespaces)
        self.header = dict()
        for element in root.find(self._get_ns('Header'), self.namespaces):
            self.header[element.tag.rpartition('}')[2]] = element.text
        self._extra_steps(root, timeout)

    def _extra_steps(self, root: ET.Element, timeout: int) -> None:
        '''Sets up subclasses. use this instead of super().__init__ because
        do not want to carry around root which could be big. Probably a better
        way.'''

    def __repr__(self) -> str:
        return f'{type(self).__name__}(\'{self.uri}\')'

    def __str__(self) -> str:
        try:
            return f'{type(self).__name__} File: {self.name}.'
        except AttributeError:
            return f'{type(self).__name__} File.'

    def __hash__(self) -> int:
        return hash(self.uri)

    def __eq__(self, other) -> bool:
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
            'Enumeration': 'http://www.sdmx.org/resources/'
                           'sdmxml/schemas/v2_1/structure',
            'Codelist': 'http://www.sdmx.org/resources/'
                        'sdmxml/schemas/v2_1/structure',
            'Code': 'http://www.sdmx.org/resources/sdmxml'
                    '/schemas/v2_1/structure',
            'Name': 'http://www.sdmx.org/resources/sdmxml/'
                    'schemas/v2_1/common',
            'DataStructureComponents': 'http://www.sdmx.org/resources/'
                                       'sdmxml/schemas/v2_1/structure'
            }
        for ns_prefix, ns_uri in self.namespaces.items():
            if ns_uri == ns_dict[element_name]:
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


class _ValidateSeriesWithDSD():
    '''Private class used to mix these functions into DSD and series.'''

    @staticmethod
    def _validate_series(dsd, series) -> bool:
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
    def _validate_series_details(dsd, series) -> Dict[str, bool]:
        '''Returns a dict with true or false for each dimension on the DSD.'''
        conformingdict = dict()
        for dim, code_list in dsd.dimensions.items():
            if dim in series.metadata.keys():
                conformingdict[dim] = series.metadata[dim] in code_list.codes()
            else:
                conformingdict[dim] = False
        return conformingdict

    @staticmethod
    def _validate_series_dimnension(dsd, series, dimension: str) -> bool:
        '''Returns true if series has a valid entery in the
        supplied dimension.
        '''
        try:
            return series.metadata[dimension] in dsd.dimensions[dimension].codes()
        except KeyError:
            raise Exception(f'Invalid dimension: {dimension}')

    @staticmethod
    def _name_from_metadata(dsd, series) -> Dict[str, str]:
        '''Returns a dictonary with names for dimensions and values.'''
        human_readable = dict()
        for dimension_id, code_list in dsd.dimensions.items():
            human_readable[dimension_id] = code_list.name_from_code(series.metadata[dimension_id])
            # Not possible to resolve dimension id for example
            # counterpart area (the dimension name) in ECOFIN
        return human_readable


class DSD(_ValidateSeriesWithDSD, SDMX):
    '''The DSD class is based on the SDMX class. It adds a dict
    of dimension ids and code lists. This class has methods which can be used
    validate a series or return the names for both dimensions and values for
    a given series.
    '''

    __slots__ = ['name', 'dimensions']

    def _extra_steps(self, root: ET.Element, timeout: int) -> None:
        self.name = root.find(f'.//{self._get_ns("DataStructure")}'
                              f'/{self._get_ns("Name")}', self.namespaces
                              ).text
        self.dimensions = dict()
        for dimension in root.find(f'.//{self._get_ns("DataStructure")}/'
                                   f'{self._get_ns("DataStructureComponents")}'
                                   f'/{self._get_ns("DimensionList")}',
                                   self.namespaces):
            for dimension_ref in dimension.findall(f'.//{self._get_ns("Enumeration")}'
                                                   f'/Ref',
                                                   self.namespaces):
                cl_url = self._generate_cl_url(dimension_ref)
                # TODO: This isn't going to work for imbeded codelists
                # such as those produced by EcOS. We still have a list
                # Name and indicator list are all that are important.
                # Maybe make a structure to hold them and use CodeList class
                # to populate if URL otherwise populate from DSD XML.
                self.dimensions[dimension.attrib['id']] = CodeList(cl_url, timeout=timeout)
        # TODO: Low priority, add attributes and measures

    def __len__(self) -> int:
        return len(self.dimensions)

    @staticmethod
    def _generate_cl_url(cl_ref: ET.Element) -> str:
        cl_agency = cl_ref.attrib['agencyID']
        cl_version = cl_ref.attrib['version']
        cl_id = cl_ref.attrib['id']
        ws_endpoint = ENTRY_URL.get(cl_agency, ENTRY_URL['SDMX'])
        return f'{ws_endpoint}codelist/{cl_agency}/{cl_id}/{cl_version}/?' \
               f'format=sdmx-2.1&detail=full&references=none'

    def validate_series(self, series) -> bool:
        '''Returns true if all dimensions for DSD are present and populated
        with valid values.
        '''
        return self._validate_series(self, series)

    def validate_series_details(self, series) -> Dict[str, bool]:
        '''Returns a dict with true or false for each dimension.'''
        return self._validate_series_details(self, series)

    def validate_series_dimnension(self, series, dimension: str) -> bool:
        '''Returns true if series has a valid entery in the given dimension.'''
        return self._validate_series_dimnension(self, series, dimension)

    def name_from_metadata(self, series) -> Dict[str, str]:
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
            self.indicators[code.attrib['id']] = code.find(f'.//{self._get_ns("Name")}',
                                                           self.namespaces).text

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
            raise Exception('Invalid code: {code}')


class SDMXTimeseries(_ValidateSeriesWithDSD):
    '''SDMXTimeseries class allows for basic manipulation for timeseries.'''

    __slots__ = ['metadata', 'frequency', 'scale', 'data', 'invalid_data',
                 'time_period_coverage']

    def __init__(self, xmlseries: ET.Element) -> None:
        self.metadata = xmlseries.attrib.copy()
        self.frequency = self.metadata.get('FREQ', 'UNKNOWN')[0]
        self.scale = self.metadata.get('UNIT_MULT', 'UNKNOWN')[0]
        self.data = dict()
        self.invalid_data = dict()
        for obs in xmlseries.iter('Obs'):
            try:
                self.data[obs.attrib['TIME_PERIOD']] = float(obs.attrib['OBS_VALUE'])
            except ValueError:
                self.invalid_data[obs.attrib['TIME_PERIOD']] = obs.attrib['OBS_VALUE']
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
        return self.time_period_coverage[-1]

    @property
    def has_invalid_data(self) -> bool:
        '''Returns last date with has data.'''
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
        self.dsd = self._generate_dsd_url(self.namespaces)
        self.series = []
        if (self.version == '2.1' or self.version == '2.0'):
            for series in root.iter('Series'):
                self.series.append(SDMXTimeseries(series))

    def __len__(self) -> int:
        return len(self.series)

    @staticmethod
    def _generate_dsd_url(namespaces: Dict[str, str]) -> str:
        phrase = 'urn:sdmx:org.sdmx.infomodel.datastructure.datastructure='
        for uri in namespaces.values():
            if uri.lower().startswith(phrase):
                agency_id = uri[len(phrase):].split(':')[0]
                ws_endpoint = ENTRY_URL.get(agency_id, ENTRY_URL['SDMX'])
                dsd = uri[len(phrase):].split(':')[1]
                dsd_id = dsd[:dsd.find('(')]
                dsd_version = dsd[dsd.find('(')+1:dsd.rfind(')')]
                return f'{ws_endpoint}datastructure/{agency_id}/{dsd_id}/' \
                       f'{dsd_version}/?format=sdmx-2.1&' \
                       f'detail=full&references=none'
        return 'UNKNOWN'

    def generate_dsd(self) -> DSD:
        '''Returns an instance of the DSD used for this SDMX file.'''
        if self.dsd == 'UNKNOWN':
            raise Exception('Could not parse DSD url.')
        else:
            return DSD(self.dsd)

    def dimensions(self) -> Tuple[str]:
        '''Returns all dimensions used in root.  All series are looped over
        just in case their is an inconsistency in the file
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
            raise Exception('Invalid dimension: {dimension}')

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
