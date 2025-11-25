# Hangman Game

A classic hangman word-guessing game with multiple categories.

## How to Play

### Starting the Game

1. Open your command line/terminal
2. Navigate to the game directory:
   ```
   cd c:\Users\eg....
   ```
3. Run the game:
   ```
   python hangman.py
   ```

### Game Rules

1. **Choose a Category**: Select from 5 categories:
   - Famous People
   - Animals
   - Car Brands
   - Professions
   - Movie Names

2. **Guess Letters**: The game will display underscores for each letter in the word/phrase. Spaces between words are visible.

3. **Make Your Guesses**: 
   - Type a single letter and press Enter
   - If the letter is in the word, it will be revealed
   - If the letter is not in the word, the hangman starts to build

4. **Win or Lose**:
   - **Win**: Guess all the letters before making 5 mistakes
   - **Lose**: Make 5 incorrect guesses and the hangman is complete

5. **Play Again**: After each game, you can choose to play again or exit

### Example Gameplay

```
Choose a category:
1. Famous People
2. Animals
3. Car Brands
4. Professions
5. Movie Names

Enter your choice (1-5): 2

Word: _ _ _ _ _ _ _ _

Guess a letter: e
✓ Correct! The letter 'e' is in the word!

Word: e _ e _ _ _ _ _
```

## Features

- 10 words/phrases per category
- Visual hangman drawing that builds with each mistake
- Tracks guessed letters to avoid duplicates
- Clear display of current word state
- 5 mistakes maximum before game over

Enjoy playing! 🎮
