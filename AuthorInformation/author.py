"""

Create an author. This is the information we have

0	Blank - this is just the id
1	ORCID
2	acronym
3	Last Name
4	First Name
5	Middle Initial
6	Email Address
7	Department
8	Institution
9	Street Address
10	City
11	State / Region
12	ZIP / Postcode
13	Country
14	Order in the list
15	Contribution
16
17	Department
18	Institution
19	Street Address
20	City
21	State / Region
22	ZIP / Postcode
23	Country





"""

import os
import sys
import argparse


class Address:
    """
    The address of an institution
    """

    def __init__(self, verbose=False):
        """
        Create a new address
        :param verbose: add some additional output
        :type verbose: bool
        """
        self.department = None
        self.institution = None
        self.street = None
        self.city = None
        self.state = None
        self.zip = None
        self.country = None
        self.verbose = verbose

    def __eq__(self, other):
        """
        Two address are equal if they have the same street, city, state, zip, and country

        :param other: The other address
        :type other: Address
        :return: If they are equal
        :rtype: bool
        """
        if isinstance(other, Address):
            return (self.street, self.city, self.state, self.zip) \
                == (other.street, other.city, other.state, other.zip)
        else:
            return NotImplemented


    def __cmp__(self, other):
        """
        Compare whether two things are the same.

        :param other: The other Address
        :type other: Address
        :return: An int, zero if they are the same
        :rtype: int
        """
        if isinstance(other, Address):
            if __eq__(other):
                return 0
            else:
                return 1
        else:
            return NotImplemented


    def __ne__(self, other):
        """
        Are these not equal?

        :param other: The other Address
        :type other: Address
        :return: If they are not equal
        :rtype: bool
        """
        result = self.__eq__(other)
        if result is NotImplemented:
            return result
        return not result


    def __hash__(self):
        """
        The hash function is based on the name of the compound.

        :rtype: int
        """
        return hash(self.__str__)


    def __str__(self):
        """
        The to string function.
        We will develop this as we go along!
        :rtype: str
        """

        toreturn = ""

        if self.department:
            toreturn = "{}, ".format(self.department)
        if self.institution:
            toreturn += "{}, ".format(self.institution)
        if self.street:
            toreturn += "{}, ".format(self.street)
        if self.city:
            toreturn += "{}, ".format(self.city)
        if self.state:
            toreturn += "{}, ".format(self.state)
        if self.zip:
            toreturn += "{}, ".format(self.zip)
        if self.country:
            toreturn += "{}".format(self.country)
        return toreturn

    def get_address(self):
        """
        Get the address
        :return:
        """
        return self.__str__()

    def is_valid(self):
        """
        Determines if this is a valid address. We need at least institution, street, city, zip and country
        :return: whether it is valid
        :rtype : bool
        """

        if self.verbose:
            if not self.institution:
                sys.stderr.write("Invalid address. No institution found\n")
            if not self.street:
                sys.stderr.write("Invalid address. No street found\n")
            if not self.city:
                sys.stderr.write("Invalid address. No city found\n")
            if not self.country:
                sys.stderr.write("Invalid address. No country found\n")


        if (self.institution and self.street and self.city and self.country):
            return True
        return False


class Author:
    """
    An author is hopefully a person
    """

    def __init__(self, abbreviation, verbose=False):
        """
        Create a new author
        :param abbreviation: a two to three letter acronym for the author. This must be unique among all authors
        :type abbreviation: str
        :param verbose: print more output
        :type verbose: bool
        """
        self.abbreviation = abbreviation
        self.orcid = None
        self.lastname = None
        self.lastnamelower = None
        self.firstname = None
        self.firstnamelower = None
        self.middleinitial = None
        self.email = None
        self.primaryaddress = Address()
        self.secondaryaddress = Address()
        self.contribution = None
        self.funding = None
        self.order = 100 # this is high by default but we overwrite it if we have a value
        self.verbose = verbose

    def __eq__(self, other):
        """
        Two authors are equal if they have the same abbreviation

        :param other: The other author
        :type other: Author
        :return: If they are equal
        :rtype: bool
        """
        if isinstance(other, Author):
            return self.abbreviation == other.abbreviation
        else:
            return NotImplemented

    def __cmp__(self, other):
        """
        Compare whether two things are the same.

        :param other: The other author
        :type other: Author
        :return: An int greater than one if we are bigger, less than one if we are smaller, zero if they are the same
        :rtype: int
        """
        if isinstance(other, Author):

            if self.order and other.order:
                if self.order < other.order:
                    return -10
                elif self.order > other.order:
                    return 10
                else:
                    return self.__str__().lower().__cmp__(other.__str__().lower())

            if self.order:
                return -10

            if other.order:
                return 10

            return self.__str__().lower().__cmp__(other.__str__().lower())
        else:
            return NotImplemented

    def __ne__(self, other):
        """
        Are these not equal?

        :param other: The other author
        :type other: Author
        :return: If they are not equal
        :rtype: bool
        """
        result = self.__eq__(other)
        if result is NotImplemented:
            return result
        return not result

    def __hash__(self):
        """
        The hash function is based on the name of the compound.

        :rtype: int
        """
        return hash(self.abbreviation)

    def __str__(self):
        """
        The to string function.
        We will develop this as we go along!
        :rtype: str
        """
        toreturn = "{}, {}".format(self.lastname, self.firstname)
        if self.middleinitial:
            toreturn = "{} {}.".format(toreturn, self.middleinitial)
        return toreturn

    def get_name(self):
        """
        get the authors name
        :return:
        """
        return self.__str__()


    def get_primary_address(self):
        """
        Return the primary address as a string
        :return:
        """
        if not self.primaryaddress.is_valid(self.verbose):
            if self.verbose:
                sys.stderr.write("Primary address not valid\n")
            return None
        else:
            return self.primaryaddress.__str__()

    def get_secondary_address(self):
        """
        Return the primary address as a string
        :return:
        """
        if not self.secondaryaddress.is_valid(self.verbose):
            if self.verbose:
                sys.stderr.write("Primary address not valid\n")
            return None
        else:
            return self.secondaryaddress.__str__()


    def is_valid(self):
        """
        valid authors have firstname, lastname, valid address, and contribution
        :return: bool
        """

        self.primaryaddress.verbose = self.verbose
        self.secondaryaddress.verbose = self.verbose

        if self.verbose:
            if not self.firstname:
                sys.stderr.write("No first name for {}\n".format(self.abbreviation))
            if not self.lastname:
                sys.stderr.write("No last name for {}\n".format(self.abbreviation))
            if not self.contribution:
                sys.stderr.write("No contribution for {}\n".format(self.abbreviation))

        if (self.firstname and self.lastname and self.contribution and self.primaryaddress.is_valid()):
            return True
        return False