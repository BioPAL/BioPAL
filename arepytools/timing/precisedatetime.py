# SPDX-FileCopyrightText: Aresys S.r.l. <info@aresys.it>
# SPDX-License-Identifier: MIT

"""
Precise Date-Time module
------------------------
"""

import datetime
import math
import copy
import re
from dateutil import parser, tz


_ISO_DATE_SEPARATOR = "-"
_ISO_TIME_SEPARATOR = ":"
_ISO_FRACTION_REGEX = re.compile('[\\.,]([0-9]+)')


def _parse_isoformat_date(date_string):
    """Parse string formatted according to ISO, returns date components and remaining chars to parse."""
    # parse year
    if len(date_string) < 4:
        raise ValueError("Incomplete year component")
    year = int(date_string[:4])
    date_string = date_string[4:]
    if not date_string:
        return (year, 1, 1), date_string

    # skip separator if present
    has_separator = date_string[0] == _ISO_DATE_SEPARATOR
    if has_separator:
        date_string = date_string[1:]

    # parse month
    if len(date_string) < 2:
        raise ValueError("Incomplete month component")
    month = int(date_string[:2])
    date_string = date_string[2:]
    if not date_string:
        if has_separator:
            return (year, month, 1), date_string
        raise ValueError("Invalid ISO date")

    # skip separator if present
    if has_separator:
        if date_string[0] != _ISO_DATE_SEPARATOR:
            raise ValueError("Unexpected date separator")
        date_string = date_string[1:]

    # parse day
    if len(date_string) < 2:
        raise ValueError("Incomplete day component")
    day = int(date_string[:2])
    date_string = date_string[2:]

    return (year, month, day), date_string


def _isoformat_date(year, month, day):
    return "{:04d}-{:02d}-{:02d}".format(year, month, day)


def _parse_isoformat_tz(time_string):
    if not time_string:
        return None, time_string

    if time_string[0] == "z" or time_string[0] == "Z":
        return tz.UTC, time_string[1:]

    if time_string[0] == "+" or time_string[0] == "-":
        raise ValueError("Unsupported timezone")

    return None, time_string


def _parse_isoformat_time(time_string):
    """Parse string formatted according to ISO, returns time components and remaining chars to parse."""
    hour, minute, second, picosecond, timezone = 0, 0, 0, 0, None

    # parse hour
    if len(time_string) < 2:
        raise ValueError("Incomplete hour component")
    hour = int(time_string[:2])
    time_string = time_string[2:]
    if not time_string:
        return (hour, minute, second, picosecond, timezone), time_string

    # time zone
    timezone, time_string = _parse_isoformat_tz(time_string)
    if timezone is not None:
        return (hour, minute, second, picosecond, timezone), time_string

    # skip separator
    has_separator = time_string[0] == _ISO_TIME_SEPARATOR
    if has_separator:
        time_string = time_string[1:]

    # parse minute
    if len(time_string) < 2:
        raise ValueError("Incomplete minute component")
    minute = int(time_string[:2])
    time_string = time_string[2:]
    if not time_string:
        return (hour, minute, second, picosecond, timezone), time_string

    # time zone
    timezone, time_string = _parse_isoformat_tz(time_string)
    if timezone is not None:
        return (hour, minute, second, picosecond, timezone), time_string

    # skip separator
    if has_separator:
        if time_string[0] != _ISO_TIME_SEPARATOR:
            raise ValueError("Unexpected time separator")
        time_string = time_string[1:]

    # parse second
    if len(time_string) < 2:
        raise ValueError("Incomplete second component")
    second = int(time_string[:2])
    time_string = time_string[2:]
    if not time_string:
        return (hour, minute, second, picosecond, timezone), time_string

    # time zone
    timezone, time_string = _parse_isoformat_tz(time_string)
    if timezone is not None:
        return (hour, minute, second, picosecond, timezone), time_string

    # parse fraction of second
    fraction_match = _ISO_FRACTION_REGEX.match(time_string)
    picosecond = 0
    if fraction_match:
        fraction = fraction_match.group(1)
        digits_count = len(fraction)

        scale_factor = 10 ** (12 - digits_count)
        picosecond = int(fraction) * scale_factor

        time_string = time_string[len(fraction_match.group()):]

    # time zone
    timezone, time_string = _parse_isoformat_tz(time_string)

    return (hour, minute, second, picosecond, timezone), time_string


def _isoformat_time(hour, minute, second, picosecond, timespec="auto"):
    fmt_specs = {
        "hours": "{:02d}",
        "minutes": "{:02d}:{:02d}",
        "seconds": "{:02d}:{:02d}:{:02d}",
        "milliseconds": "{:02d}:{:02d}:{:02d}.{:03d}",
        "microseconds": "{:02d}:{:02d}:{:02d}.{:06d}",
        "nanoseconds": "{:02d}:{:02d}:{:02d}.{:09d}",
        "picoseconds": "{:02d}:{:02d}:{:02d}.{:012d}",
    }

    fraction_of_second = int(picosecond)
    if timespec == "auto":
        timespec = "picoseconds" if picosecond else "seconds"
    elif timespec == "milliseconds":
        fraction_of_second //= 1_000_000_000
    elif timespec == "microseconds":
        fraction_of_second //= 1_000_000
    elif timespec == "nanoseconds":
        fraction_of_second //= 1_000

    try:
        fmt = fmt_specs[timespec]
    except KeyError as exc:
        raise ValueError("Unknown timespec value") from exc

    return fmt.format(hour, minute, second, fraction_of_second)


class PreciseDateTime:
    """Precise Date-Time class

    Precision: 1e-12 (picoseconds)

    Standard format: "DD-MMM-YYYY hh:mm:ss.pppppppppppp"

    * DD is the day
    * MMM is the month string (JAN, FEB, MAR, APR, MAY, JUN, JUL, AUG, SEP, OCT, NOV, DEC)
    * YYYY is the year
    * hh are the hours (00->23)
    * mm are the minutes (00->59)
    * ss are the seconds (00->59)
    * pppppppppppp are the picoseconds

    Standard reference date: 01/01/1985 00:00:00"""

    _MONTH_ABBREVIATED_NAME_DIRECTIVE = '%b'
    _MONTH_ABBREVIATED_NAMES = (
        'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC')
    _STRING_FORMAT = '%d-%b-%Y %H:%M:%S.'
    _REFERENCE_DATETIME = datetime.datetime(year=1985, month=1, day=1)
    _TIME_DIFF_REFERENCE_FROM_1985 = _REFERENCE_DATETIME - datetime.datetime(year=1985, month=1, day=1)
    _PRECISION = 1e-12  # Precision of the decimal part
    _SECONDS_IN_A_DAY = 24 * 60 * 60

    def __init__(self, seconds=0., picoseconds=0.):
        """Initialize the object at the time point defined by adding the specified input seconds and picoseconds to the
        reference date.

        :param seconds: number of seconds from the reference date
        :param picoseconds: number of picoseconds to add to the specified seconds
        """
        if seconds < 0 or picoseconds < 0:
            raise ValueError("Seconds and picoseconds must be non-negative")

        self._set_state(seconds, picoseconds)

    def _set_state(self, seconds, picoseconds):
        """Set the object at the time point defined by adding the specified input seconds and picoseconds to the
        reference date.

        :param seconds: number of seconds from the reference date
        :param picoseconds: number of picoseconds to add to the specified seconds
        """
        seconds_fraction = seconds - int(seconds)
        picoseconds_in_seconds_fraction = seconds_fraction / self._PRECISION
        tot_picoseconds = picoseconds + picoseconds_in_seconds_fraction

        seconds_adj = math.floor(tot_picoseconds * self._PRECISION)
        normalized_seconds = int(seconds) + int(seconds_adj)
        normalized_picoseconds = float(tot_picoseconds) % (1 / self._PRECISION)

        if normalized_seconds < 0:
            raise ValueError("The specified time is before the reference date")

        assert normalized_seconds >= 0 and 0 <= normalized_picoseconds < 1 / self._PRECISION

        self._seconds = normalized_seconds
        self._picoseconds = normalized_picoseconds

        assert isinstance(self._seconds, int)

    def __iadd__(self, seconds):
        """Add the input seconds to the current time point.

        :param seconds: number of seconds to add to the current time point

        :return: self
        """
        seconds_fraction = seconds - int(seconds)
        self._set_state(self._seconds + int(seconds), self._picoseconds + seconds_fraction / self._PRECISION)
        return self

    def __isub__(self, seconds):
        """Subtract the input seconds from the current time point.

        :param seconds: number of seconds to subtract from the current time point

        :return: self
        """
        seconds_fraction = seconds - int(seconds)
        self._set_state(self._seconds - int(seconds), self._picoseconds - seconds_fraction / self._PRECISION)
        return self

    def __add__(self, seconds):
        """Return the sum between the current time point and the specified input seconds.

        :param seconds: number of seconds to add to the current time point

        :return: a new PreciseDateTime object initialized to the resulting time point
        """
        ret_obj = copy.deepcopy(self)
        ret_obj += seconds
        return ret_obj

    __radd__ = __add__

    def __sub__(self, other):
        """Return the difference between the current time point and the specified input parameter (seconds or another
        PreciseDateTime object).

        :param other: number of seconds or PreciseDateTime object to subtract from the current time point

        :return: if the input parameter is a PreciseDateTime object, the difference in seconds between the two time
        points; otherwise, a new PreciseDateTime object initialized to the resulting time point.
        """
        if isinstance(other, PreciseDateTime):
            seconds_fraction = (self._picoseconds - other._picoseconds) * self._PRECISION
            return self._seconds - other._seconds + seconds_fraction
        else:
            ret_obj = copy.deepcopy(self)
            ret_obj -= other
            return ret_obj

    def __repr__(self):
        assert isinstance(self._seconds, int)
        assert 0 <= self._picoseconds < 1 / self._PRECISION
        absolute_datetime = self._REFERENCE_DATETIME + datetime.timedelta(0, self._seconds)
        tmp_str = '{:0>12d}'.format(int(self._picoseconds))

        # Replacing month abbreviated name directive with english abbreviated month name
        # to be locale independent
        month_id = int(absolute_datetime.strftime('%m')) - 1
        month_name = self._MONTH_ABBREVIATED_NAMES[month_id]
        updated_string_format = self._STRING_FORMAT.replace(self._MONTH_ABBREVIATED_NAME_DIRECTIVE,
                                                            month_name)

        return absolute_datetime.strftime(updated_string_format) + tmp_str

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (other._seconds == self._seconds) and (other._picoseconds == self._picoseconds)
        else:
            return NotImplemented

    def __ne__(self, other):
        if isinstance(other, self.__class__):
            return not self == other
        else:
            return NotImplemented

    def __lt__(self, other):
        if isinstance(other, self.__class__):
            return self._seconds < other._seconds or (
                self._seconds == other._seconds and self._picoseconds < other._picoseconds)
        else:
            return NotImplemented

    def __le__(self, other):
        if isinstance(other, self.__class__):
            return self < other or self == other
        else:
            return NotImplemented

    def __gt__(self, other):
        if isinstance(other, self.__class__):
            return not self <= other
        else:
            return NotImplemented

    def __ge__(self, other):
        if isinstance(other, self.__class__):
            return not self < other
        else:
            return NotImplemented

    @classmethod
    def get_precision(cls):
        """Date-time representation precision.

        :return: precision of the date-time representation in seconds
        """
        return cls._PRECISION

    @classmethod
    def get_reference_datetime(cls):
        """Reference date-time.

        :return: reference date-time as a datetime object
        """
        return cls._REFERENCE_DATETIME

    @property
    def year(self):
        """Year associated to the current time point."""
        assert 0 <= self._picoseconds < 1 / self._PRECISION
        absolute_datetime = self._REFERENCE_DATETIME + datetime.timedelta(0, self._seconds)
        return int(absolute_datetime.strftime('%Y'))

    @property
    def month(self):
        """Month associated to the current time point."""
        assert 0 <= self._picoseconds < 1 / self._PRECISION
        absolute_datetime = self._REFERENCE_DATETIME + datetime.timedelta(0, self._seconds)
        return int(absolute_datetime.strftime('%m'))

    @property
    def day_of_the_month(self):
        """Day of the month associated to the current time point."""
        assert 0 <= self._picoseconds < 1 / self._PRECISION
        absolute_datetime = self._REFERENCE_DATETIME + datetime.timedelta(0, self._seconds)
        return int(absolute_datetime.strftime('%d'))

    @property
    def hour_of_day(self):
        """Hour of the day associated to the current time point."""
        assert 0 <= self._picoseconds < 1 / self._PRECISION
        absolute_datetime = self._REFERENCE_DATETIME + datetime.timedelta(0, self._seconds)
        return int(absolute_datetime.strftime('%H'))

    @property
    def minute_of_hour(self):
        """Minute of the hour associated to the current time point."""
        assert 0 <= self._picoseconds < 1 / self._PRECISION
        absolute_datetime = self._REFERENCE_DATETIME + datetime.timedelta(0, self._seconds)
        return int(absolute_datetime.strftime('%M'))

    @property
    def second_of_minute(self):
        """Second of the minute associated to the current time point."""
        assert 0 <= self._picoseconds < 1 / self._PRECISION
        absolute_datetime = self._REFERENCE_DATETIME + datetime.timedelta(0, self._seconds)
        return int(absolute_datetime.strftime('%S'))

    @property
    def picosecond_of_second(self):
        """Picosecond of the second associated to the current time point."""
        assert isinstance(self._seconds, int)
        assert 0 <= self._picoseconds < 1 / self._PRECISION
        return self._picoseconds

    @property
    def fraction_of_day(self):
        """Fraction of the day associated to the current time point."""
        assert 0 <= self._picoseconds < 1 / self._PRECISION
        seconds_from_day_start = self._seconds % self._SECONDS_IN_A_DAY
        return (seconds_from_day_start + self._picoseconds * self._PRECISION) / self._SECONDS_IN_A_DAY

    @property
    def day_of_the_year(self):
        """Day from the first day of the year associated to the current time point."""
        assert 0 <= self._picoseconds < 1 / self._PRECISION
        absolute_datetime_first_day_of_year = datetime.datetime(self.year, 1, 1)
        absolute_datetime = self._REFERENCE_DATETIME + datetime.timedelta(0, self._seconds)
        time_diff_from_first_day_of_year = absolute_datetime - absolute_datetime_first_day_of_year
        return int(time_diff_from_first_day_of_year.days + 1)

    @property
    def sec85(self):
        """Time distance in seconds from 01/01/1985 00:00:00 to the current time point."""
        return self._TIME_DIFF_REFERENCE_FROM_1985.total_seconds() + self._seconds + self._picoseconds * self._PRECISION

    def set_now(self):
        """Set the object to the current time (local timezone).

        :return: self
        """
        absolute_datetime_now = datetime.datetime.now()
        time_diff_from_reference_date = absolute_datetime_now - self._REFERENCE_DATETIME
        seconds = time_diff_from_reference_date.total_seconds()
        self._set_state(seconds, 0)
        return self

    def set_from_sec85(self, seconds):
        """Set the object to the time point defined by adding the specified input seconds to 01/01/1985 00:00:00.

        :param seconds: number of seconds from 01/01/1985 00:00:00.0

        :return: self
        """
        self._set_state(seconds - self._TIME_DIFF_REFERENCE_FROM_1985.total_seconds(), 0)
        return self

    def set_from_utc_string(self, utc_str):
        """Set the object to the time point specified by the input UTC string.

        :param utc_str: UTC string

        :return: self
        """
        try:
            seconds_fraction = float('0.' + utc_str.split('.')[1].strip())
            absolute_datetime = parser.parse(utc_str).replace(microsecond=0)
        except Exception as exc:
            raise InvalidUtcString(utc_str) from exc

        time_diff_from_reference_date = absolute_datetime - self._REFERENCE_DATETIME

        seconds = time_diff_from_reference_date.total_seconds()
        picoseconds = seconds_fraction / self._PRECISION

        self._set_state(seconds, picoseconds)
        return self

    def set_from_numeric_datetime(self, year, month=1, day=1, hours=0, minutes=0, seconds=0, picoseconds=0.):
        """Set the object to the time point specified by the input date and time parameters.

        :param year: year
        :param month: month, from 1 to 12
        :param day: day, from 1 to 28-31 (depending on month)
        :param hours: hours, from 0 to 23
        :param minutes: minutes, from 0 to 59
        :param seconds: seconds, from 0 to 59
        :param picoseconds: picoseconds, non-negative and less than 1e12

        :return: self
        """
        absolute_datetime = datetime.datetime(year, month, day, hours, minutes, seconds)
        if not 0 <= picoseconds < 1 / self._PRECISION:
            raise ValueError("Picoseconds must be non-negative and less than {}".format(1 / self._PRECISION))

        time_diff_from_reference_date = absolute_datetime - self._REFERENCE_DATETIME
        seconds = time_diff_from_reference_date.total_seconds()

        self._set_state(seconds, picoseconds)
        return self

    def set_from_isoformat(self, datetime_string, sep="T"):
        """Set the object to the time specified by the input ISO string.

        :param datetime_string: time specified as ISO string
        :param sep: the separator between date and time, default 'T'
        :return: self
        """
        if not isinstance(datetime_string, str):
            raise ValueError("fromisoformat: argument must be a str")

        date_string = datetime_string
        try:
            date, time_string = _parse_isoformat_date(date_string)
        except ValueError as exc:
            raise ValueError("Invalid isoformat string: {}".format(datetime_string)) from exc

        time = ()
        if time_string:
            if time_string[0] == sep:
                time_string = time_string[1:]

            try:
                time, unparsed_string = _parse_isoformat_time(time_string)

                timezone = time[-1]
                if timezone is not None and timezone != tz.UTC:
                    raise ValueError("Unsupported timezone: {}".format(timezone))
                time = time[:-1]
            except ValueError as exc:
                raise ValueError("Invalid isoformat string: {}".format(datetime_string)) from exc

            if unparsed_string:
                raise ValueError("Invalid isoformat string: {}".format(datetime_string))

        return self.set_from_numeric_datetime(*(date + time))

    @staticmethod
    def fromisoformat(datetime_string, sep="T"):
        """Construct a PreciseDateTime from a string formatted according to ISO.

        :param datetime_string: time specified as ISO string
        :param sep: the separator between date and time, default 'T'
        :return: a PreciseDateTime corresponding to the provided string
        """
        return PreciseDateTime().set_from_isoformat(datetime_string, sep=sep)

    def isoformat(self, sep="T", timespec="auto"):
        """Return the time formatted according to ISO.

        The optional argument timespec specifies the number of additional terms of the time to include. Valid options
        are 'auto', 'hours', 'minutes', 'seconds', 'milliseconds', 'microseconds', 'nanoseconds' and 'picoseconds'.

        :param sep: the separator between date and time, default 'T'
        :param timespec: the number of additional terms of the time to include
        :return: time formatted according to ISO
        """

        absolute_datetime = self._REFERENCE_DATETIME + datetime.timedelta(seconds=self._seconds)
        date = _isoformat_date(absolute_datetime.year, absolute_datetime.month, absolute_datetime.day)
        time = _isoformat_time(absolute_datetime.hour, absolute_datetime.minute, absolute_datetime.second,
                               self._picoseconds, timespec=timespec)
        return "{}{:s}{}Z".format(date, sep, time)


class InvalidUtcString(ValueError):
    """
    Invalid UTC string exception class
    """

    def __init__(self, value):
        super().__init__("Invalid UTC string: '{}'".format(value))
