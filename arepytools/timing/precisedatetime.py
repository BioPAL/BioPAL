"""
Precise Date-Time module
------------------------
"""

import datetime
import math
import copy
from dateutil import parser


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
        return (absolute_datetime.strftime(self._STRING_FORMAT) + tmp_str).upper()

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
        except Exception as e:
            raise InvalidUtcString(utc_str, e)

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


class InvalidUtcString(ValueError):
    """
    Invalid UTC string exception class
    """

    def __init__(self, value, original_exception):
        super().__init__("Invalid UTC string: '{}'".format(value))
        self.original_exception = original_exception
