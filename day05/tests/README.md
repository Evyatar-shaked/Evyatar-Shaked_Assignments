# Tests for Hangman Game

This folder contains comprehensive unit tests for the `hangman.py` game.

## Test Coverage

The test suite includes:

### 1. **TestHangmanWordBanks**
- Verifies all expected categories exist
- Ensures each category has sufficient words
- Validates word bank structure

### 2. **TestHangmanStages**
- Tests the visual hangman stages count
- Validates stage content

### 3. **TestGetRandomWord**
- Tests random word selection functionality
- Ensures words are from correct categories
- Verifies use of random.choice

### 4. **TestDisplayWordState**
- Tests word display with no letters guessed
- Tests partial letter revelation
- Tests complete word revelation
- Tests multi-word phrases with spaces
- Tests case-insensitive matching

### 5. **TestGetCategoryChoice**
- Tests valid category selections (1-5)
- Tests invalid input handling
- Tests input validation loop

### 6. **TestDisplayCategoryMenu**
- Tests menu display functionality
- Validates all categories are shown

### 7. **TestPlayGameIntegration**
- Integration test for winning scenario
- Integration test for losing scenario
- Tests duplicate guess handling

### 8. **TestMainFunction**
- Tests single game execution
- Tests multiple game loop
- Tests exit functionality

## Running the Tests

From the `day05` directory, run:

```bash
python -m pytest tests/test_hangman.py -v
```

Or using unittest:

```bash
python -m unittest tests.test_hangman
```

## Test Statistics

- **Total Test Classes**: 8
- **Total Test Cases**: 30+
- **Coverage**: All major functions and game logic paths

## Notes

- Tests use mocking for user input and random selections
- Integration tests simulate complete game scenarios
- All tests are independent and can run in any order
