import itertools

from order2_debruijn_arithmetic import *

i_range = 36
sequence = "00110212203132330414243440515253545506162636465660717273747576770818283848586878809192939495969798990A1A2A3A4A5A6A7A8A9AA0B1B2B3B4B5B6B7B8B9BABB0C1C2C3C4C5C6C7C8C9CACBCC0D1D2D3D4D5D6D7D8D9DADBDCDD0E1E2E3E4E5E6E7E8E9EAEBECEDEE0F1F2F3F4F5F6F7F8F9FAFBFCFDFEFF0G1G2G3G4G5G6G7G8G9GAGBGCGDGEGFGG0H1H2H3H4H5H6H7H8H9HAHBHCHDHEHFHGHH0I1I2I3I4I5I6I7I8I9IAIBICIDIEIFIGIHII0J1J2J3J4J5J6J7J8J9JAJBJCJDJEJFJGJHJIJJ0K1K2K3K4K5K6K7K8K9KAKBKCKDKEKFKGKHKIKJKK0L1L2L3L4L5L6L7L8L9LALBLCLDLELFLGLHLILJLKLL0M1M2M3M4M5M6M7M8M9MAMBMCMDMEMFMGMHMIMJMKMLMM0N1N2N3N4N5N6N7N8N9NANBNCNDNENFNGNHNINJNKNLNMNN0O1O2O3O4O5O6O7O8O9OAOBOCODOEOFOGOHOIOJOKOLOMONOO0P1P2P3P4P5P6P7P8P9PAPBPCPDPEPFPGPHPIPJPKPLPMPNPOPP0Q1Q2Q3Q4Q5Q6Q7Q8Q9QAQBQCQDQEQFQGQHQIQJQKQLQMQNQOQPQQ0R1R2R3R4R5R6R7R8R9RARBRCRDRERFRGRHRIRJRKRLRMRNRORPRQRR0S1S2S3S4S5S6S7S8S9SASBSCSDSESFSGSHSISJSKSLSMSNSOSPSQSRSS0T1T2T3T4T5T6T7T8T9TATBTCTDTETFTGTHTITJTKTLTMTNTOTPTQTRTSTT0U1U2U3U4U5U6U7U8U9UAUBUCUDUEUFUGUHUIUJUKULUMUNUOUPUQURUSUTUU0V1V2V3V4V5V6V7V8V9VAVBVCVDVEVFVGVHVIVJVKVLVMVNVOVPVQVRVSVTVUVV0W1W2W3W4W5W6W7W8W9WAWBWCWDWEWFWGWHWIWJWKWLWMWNWOWPWQWRWSWTWUWVWW0X1X2X3X4X5X6X7X8X9XAXBXCXDXEXFXGXHXIXJXKXLXMXNXOXPXQXRXSXTXUXVXWXX0Y1Y2Y3Y4Y5Y6Y7Y8Y9YAYBYCYDYEYFYGYHYIYJYKYLYMYNYOYPYQYRYSYTYUYVYWYXYY0Z1Z2Z3Z4Z5Z6Z7Z8Z9ZAZBZCZDZEZFZGZHZIZJZKZLZMZNZOZPZQZRZSZTZUZVZWZXZYZZ0"


class TestSanity:
    def test_rho2(self):
        assert rho2((2, 4)) == 20
        assert rho2((3, 1)) == 10

    def test_inverse_rho2(self):
        assert inverse_rho2(20) == (2, 4)
        assert inverse_rho2(10) == (3, 1)

    def test_lambda2(self):
        assert lambda2((2, 4)) == (4, 2, 0)
        assert lambda2((3, 1)) == (3, 1, 1)
    
    def test_inverse_lambda2(self):
        assert inverse_lambda2(2, 1, 0) == (1, 2)
        assert inverse_lambda2(2, 2, 0) == (2, 0)
    
    def test_carry_down(self):
        assert carry_down(7, -9, 1) == (5, 3, 1)

    def test_addition(self):
        assert add_states2((2, 4), (3, 1)) == (5, 3)

    def test_carry_up(self):
        assert carry_up(4, 16, 0) == (6, 6, 0)
    
    def test_multiplication(self):
        assert mul_states2((1, 2), (2, 0)) == (6, 0)


class TestEdgeCases:
    def test_null_state(self):
        assert rho2((0, 0)) == 0
        assert inverse_rho2(0) == (0, 0)
    
    def test_rho2_terminal_states(self):
        assert rho2((11, 11)) == 142
        assert rho2((11, 0)) == 143

    def test_inverse_rho2_terminal_states(self):
        assert inverse_rho2(62) == (7, 7)
        assert inverse_rho2(63) == (7, 0)
    
    def test_lambda2_terminal_states(self):
        assert lambda2((5, 5)) == (5, 5, 1)
        assert lambda2((5, 0)) == (5, 5, 0)

    def test_inverse_lambda2_terminal_states(self):
        assert inverse_lambda2(3, 3, 1) == (3, 3)
        assert inverse_lambda2(3, 3, 0) == (3, 0)


class TestAll:
    def test_rho2(self):
        for a, b in itertools.product(range(i_range), range(i_range)):
            index = rho2((a, b))
            assert tuple(int(i, 36) for i in sequence[index:index + 2]) == (a, b)

    def test_inverse_rho2(self):
        for i in range(i_range ** 2):
            assert tuple(int(i, 36) for i in sequence[i:i + 2]) == inverse_rho2(i)

    def test_addition(self):
        for a, b, c, d in itertools.product(range(i_range), range(i_range), range(i_range), range(i_range)):
            x = (a, b)
            y = (c, d)
            assert add_states2(x, y) == inverse_rho2(rho2(x) + rho2(y))

    def test_multiplication(self):
        for a, b, c, d in itertools.product(range(i_range), range(i_range), range(i_range), range(i_range)):
            x = (a, b)
            y = (c, d)
            assert mul_states2(x, y) == inverse_rho2(rho2(x) * rho2(y))
