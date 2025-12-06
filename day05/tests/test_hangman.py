import unittest
from unittest.mock import patch, MagicMock
import sys
import os

# Add parent directory to path to import hangman module
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import hangman


class TestHangmanWordBanks(unittest.TestCase):
    """Test the word banks and their structure."""
    
    def test_word_banks_exist(self):
        """Test that all expected categories exist in WORD_BANKS."""
        expected_categories = ["famous people", "animals", "car brands", "professions", "movie names"]
        for category in expected_categories:
            self.assertIn(category, hangman.WORD_BANKS)
    
    def test_word_banks_not_empty(self):
        """Test that each category has words."""
        for category, words in hangman.WORD_BANKS.items():
            self.assertGreater(len(words), 0, f"Category '{category}' should not be empty")
    
    def test_word_banks_have_enough_words(self):
        """Test that each category has at least 5 words for variety."""
        for category, words in hangman.WORD_BANKS.items():
            self.assertGreaterEqual(len(words), 5, f"Category '{category}' should have at least 5 words")


class TestHangmanStages(unittest.TestCase):
    """Test the hangman visual stages."""
    
    def test_hangman_stages_count(self):
        """Test that there are exactly 6 stages (0-5 mistakes)."""
        self.assertEqual(len(hangman.HANGMAN_STAGES), 6)
    
    def test_hangman_stages_not_empty(self):
        """Test that each stage has content."""
        for stage in hangman.HANGMAN_STAGES:
            self.assertIsInstance(stage, str)
            self.assertGreater(len(stage), 0)


class TestGetRandomWord(unittest.TestCase):
    """Test the get_random_word function."""
    
    def test_get_random_word_returns_string(self):
        """Test that get_random_word returns a string."""
        for category in hangman.WORD_BANKS.keys():
            word = hangman.get_random_word(category)
            self.assertIsInstance(word, str)
    
    def test_get_random_word_from_correct_category(self):
        """Test that the returned word is from the selected category."""
        for category in hangman.WORD_BANKS.keys():
            word = hangman.get_random_word(category)
            self.assertIn(word, hangman.WORD_BANKS[category])
    
    @patch('hangman.random.choice')
    def test_get_random_word_uses_random_choice(self, mock_choice):
        """Test that get_random_word uses random.choice."""
        mock_choice.return_value = "test"
        category = "animals"
        hangman.get_random_word(category)
        mock_choice.assert_called_once_with(hangman.WORD_BANKS[category])


class TestDisplayWordState(unittest.TestCase):
    """Test the display_word_state function."""
    
    def test_display_word_state_all_hidden(self):
        """Test word display when no letters are guessed."""
        word = "elephant"
        guessed_letters = set()
        result = hangman.display_word_state(word, guessed_letters)
        self.assertEqual(result, "_ _ _ _ _ _ _ _")
    
    def test_display_word_state_some_guessed(self):
        """Test word display with some letters guessed."""
        word = "elephant"
        guessed_letters = {'e', 'a'}
        result = hangman.display_word_state(word, guessed_letters)
        self.assertEqual(result, "e _ e _ _ a _ _")
    
    def test_display_word_state_all_guessed(self):
        """Test word display when all letters are guessed."""
        word = "cat"
        guessed_letters = {'c', 'a', 't'}
        result = hangman.display_word_state(word, guessed_letters)
        self.assertEqual(result, "c a t")
    
    def test_display_word_state_with_spaces(self):
        """Test word display with spaces (multi-word phrases)."""
        word = "Hello World"
        guessed_letters = {'h', 'o'}
        result = hangman.display_word_state(word, guessed_letters)
        # Spaces should be preserved with double spacing
        self.assertIn("  ", result)  # Double space for word separation
    
    def test_display_word_state_case_insensitive(self):
        """Test that guessing works case-insensitively."""
        word = "Einstein"
        guessed_letters = {'e', 'i', 'n'}
        result = hangman.display_word_state(word, guessed_letters)
        # Should reveal both 'E' and 'e', both 'i's, and both 'n's
        self.assertIn("E", result)
        self.assertIn("i", result)
        self.assertIn("n", result)


class TestGetCategoryChoice(unittest.TestCase):
    """Test the get_category_choice function."""
    
    @patch('builtins.input', return_value='1')
    def test_get_category_choice_valid_input_1(self, mock_input):
        """Test category choice with valid input '1'."""
        result = hangman.get_category_choice()
        self.assertEqual(result, "famous people")
    
    @patch('builtins.input', return_value='2')
    def test_get_category_choice_valid_input_2(self, mock_input):
        """Test category choice with valid input '2'."""
        result = hangman.get_category_choice()
        self.assertEqual(result, "animals")
    
    @patch('builtins.input', return_value='3')
    def test_get_category_choice_valid_input_3(self, mock_input):
        """Test category choice with valid input '3'."""
        result = hangman.get_category_choice()
        self.assertEqual(result, "car brands")
    
    @patch('builtins.input', return_value='4')
    def test_get_category_choice_valid_input_4(self, mock_input):
        """Test category choice with valid input '4'."""
        result = hangman.get_category_choice()
        self.assertEqual(result, "professions")
    
    @patch('builtins.input', return_value='5')
    def test_get_category_choice_valid_input_5(self, mock_input):
        """Test category choice with valid input '5'."""
        result = hangman.get_category_choice()
        self.assertEqual(result, "movie names")
    
    @patch('builtins.input', side_effect=['invalid', '0', '6', '3'])
    @patch('builtins.print')
    def test_get_category_choice_invalid_then_valid(self, mock_print, mock_input):
        """Test category choice with invalid inputs followed by valid input."""
        result = hangman.get_category_choice()
        self.assertEqual(result, "car brands")
        # Should have printed error messages for invalid inputs
        self.assertGreater(mock_print.call_count, 0)


class TestDisplayCategoryMenu(unittest.TestCase):
    """Test the display_category_menu function."""
    
    @patch('builtins.print')
    def test_display_category_menu_prints(self, mock_print):
        """Test that display_category_menu prints output."""
        hangman.display_category_menu()
        self.assertGreater(mock_print.call_count, 0)
    
    @patch('builtins.print')
    def test_display_category_menu_contains_categories(self, mock_print):
        """Test that menu contains all category options."""
        hangman.display_category_menu()
        # Combine all printed output
        all_output = ' '.join([str(call[0][0]) for call in mock_print.call_args_list])
        self.assertIn("Famous People", all_output)
        self.assertIn("Animals", all_output)
        self.assertIn("Car Brands", all_output)
        self.assertIn("Professions", all_output)
        self.assertIn("Movie Names", all_output)


class TestPlayGameIntegration(unittest.TestCase):
    """Integration tests for the play_game function."""
    
    @patch('builtins.input', side_effect=['1', 't', 'e', 's'])
    @patch('builtins.print')
    @patch('hangman.get_random_word', return_value='test')
    def test_play_game_win_scenario(self, mock_get_word, mock_print, mock_input):
        """Test a winning game scenario."""
        hangman.play_game()
        # Check that congratulations message was printed
        all_output = ' '.join([str(call[0][0]) if call[0] else '' for call in mock_print.call_args_list])
        self.assertIn("CONGRATULATIONS", all_output.upper())
    
    @patch('builtins.input', side_effect=['2', 'x', 'z', 'q', 'j', 'k', 'w'])
    @patch('builtins.print')
    @patch('hangman.get_random_word', return_value='cat')
    def test_play_game_lose_scenario(self, mock_get_word, mock_print, mock_input):
        """Test a losing game scenario."""
        hangman.play_game()
        # Check that game over message was printed
        all_output = ' '.join([str(call[0][0]) if call[0] else '' for call in mock_print.call_args_list])
        self.assertIn("LOST", all_output.upper())
    
    @patch('builtins.input', side_effect=['3', 'a', 'a', 'e'])
    @patch('builtins.print')
    @patch('hangman.get_random_word', return_value='ae')
    def test_play_game_duplicate_guess(self, mock_get_word, mock_print, mock_input):
        """Test that duplicate guesses are handled."""
        hangman.play_game()
        # Check that duplicate message was printed
        all_output = ' '.join([str(call[0][0]) if call[0] else '' for call in mock_print.call_args_list])
        self.assertIn("already guessed", all_output.lower())


class TestMainFunction(unittest.TestCase):
    """Test the main function."""
    
    @patch('builtins.input', side_effect=['1', 't', 'e', 's', 't', 'no'])
    @patch('builtins.print')
    @patch('hangman.get_random_word', return_value='test')
    def test_main_single_game(self, mock_get_word, mock_print, mock_input):
        """Test main function with one game and exit."""
        hangman.main()
        # Check that goodbye message was printed
        all_output = ' '.join([str(call[0][0]) if call[0] else '' for call in mock_print.call_args_list])
        self.assertIn("Goodbye", all_output)
    
    @patch('builtins.input', side_effect=['2', 'c', 'a', 't', 'yes', '3', 'd', 'o', 'g', 'n'])
    @patch('builtins.print')
    @patch('hangman.get_random_word', side_effect=['cat', 'dog'])
    def test_main_multiple_games(self, mock_get_word, mock_print, mock_input):
        """Test main function with multiple games."""
        hangman.main()
        # Should have called get_random_word twice
        self.assertEqual(mock_get_word.call_count, 2)


if __name__ == '__main__':
    unittest.main()
