import random

# Word banks for each category
WORD_BANKS = {
    "famous people": [
        "Albert Einstein",
        "Leonardo DiCaprio",
        "Oprah Winfrey",
        "Michael Jordan",
        "Taylor Swift",
        "Elvis Presley",
        "Marie Curie",
        "Nelson Mandela",
        "Steve Jobs",
        "Barack Obama"
    ],
    "animals": [
        "elephant",
        "giraffe",
        "penguin",
        "kangaroo",
        "butterfly",
        "crocodile",
        "dolphin",
        "cheetah",
        "octopus",
        "rhinoceros"
    ],
    "car brands": [
        "Toyota",
        "Mercedes Benz",
        "Ferrari",
        "Volkswagen",
        "Chevrolet",
        "Lamborghini",
        "Porsche",
        "Aston Martin",
        "Land Rover",
        "Rolls Royce"
    ],
    "professions": [
        "doctor",
        "engineer",
        "teacher",
        "firefighter",
        "pilot",
        "chef",
        "architect",
        "electrician",
        "veterinarian",
        "pharmacist"
    ],
    "movie names": [
        "The Godfather",
        "Forrest Gump",
        "The Matrix",
        "Pulp Fiction",
        "Jurassic Park",
        "The Lion King",
        "Star Wars",
        "Avatar",
        "The Avengers",
        "Toy Story"
    ]
}

# Hangman stages (5 mistakes = game over)
HANGMAN_STAGES = [
    """
       ------
       |    |
       |
       |
       |
       |
    --------
    """,
    """
       ------
       |    |
       |    O
       |
       |
       |
    --------
    """,
    """
       ------
       |    |
       |    O
       |    |
       |
       |
    --------
    """,
    """
       ------
       |    |
       |    O
       |   /|
       |
       |
    --------
    """,
    """
       ------
       |    |
       |    O
       |   /|\\
       |
       |
    --------
    """,
    """
       ------
       |    |
       |    O
       |   /|\\
       |   / \\
       |
    --------
    GAME OVER!
    """
]

def display_category_menu():
    """Display the category selection menu."""
    print("\n" + "="*50)
    print("WELCOME TO HANGMAN GAME!")
    print("="*50)
    print("\nChoose a category:")
    print("1. Famous People")
    print("2. Animals")
    print("3. Car Brands")
    print("4. Professions")
    print("5. Movie Names")
    print("="*50)

def get_category_choice():
    """Get and validate the user's category choice."""
    category_map = {
        "1": "famous people",
        "2": "animals",
        "3": "car brands",
        "4": "professions",
        "5": "movie names"
    }
    
    while True:
        choice = input("\nEnter your choice (1-5): ").strip()
        if choice in category_map:
            return category_map[choice]
        else:
            print("Invalid choice! Please enter a number between 1 and 5.")

def get_random_word(category):
    """Get a random word from the selected category."""
    return random.choice(WORD_BANKS[category])

def display_word_state(word, guessed_letters):
    """Display the current state of the word with guessed letters and underscores."""
    display = ""
    for char in word:
        if char == " ":
            display += "  "  # Double space for word separation visibility
        elif char.lower() in guessed_letters or not char.isalpha():
            display += char
        else:
            display += "_"
        display += " "
    return display.strip()

def play_game():
    """Main game logic."""
    # Category selection
    display_category_menu()
    category = get_category_choice()
    
    # Get random word from category
    word_to_guess = get_random_word(category)
    
    print(f"\nGreat! You chose: {category.upper()}")
    print(f"The word/phrase has {len(word_to_guess)} characters.")
    print("\nLet's start the game!")
    
    # Game variables
    guessed_letters = set()
    mistakes = 0
    max_mistakes = 5
    
    # Game loop
    while mistakes < max_mistakes:
        print("\n" + "="*50)
        print(HANGMAN_STAGES[mistakes])
        print("="*50)
        
        # Display current word state
        current_display = display_word_state(word_to_guess, guessed_letters)
        print(f"\nWord: {current_display}")
        print(f"\nGuessed letters: {', '.join(sorted(guessed_letters)) if guessed_letters else 'None'}")
        print(f"Mistakes: {mistakes}/{max_mistakes}")
        
        # Check if word is completely guessed
        if all(char.lower() in guessed_letters or not char.isalpha() for char in word_to_guess):
            print("\n" + "="*50)
            print("🎉 CONGRATULATIONS! YOU WON! 🎉")
            print("="*50)
            print(f"The word was: {word_to_guess}")
            break
        
        # Get user guess
        guess = input("\nGuess a letter: ").strip().lower()
        
        # Validate input
        if len(guess) != 1 or not guess.isalpha():
            print("Please enter a single letter!")
            continue
        
        if guess in guessed_letters:
            print(f"You already guessed '{guess}'! Try a different letter.")
            continue
        
        # Add to guessed letters
        guessed_letters.add(guess)
        
        # Check if letter is in the word
        if guess in word_to_guess.lower():
            print(f"✓ Correct! The letter '{guess}' is in the word!")
        else:
            print(f"✗ Wrong! The letter '{guess}' is not in the word.")
            mistakes += 1
    
    # Game over - player lost
    if mistakes == max_mistakes:
        print("\n" + "="*50)
        print(HANGMAN_STAGES[mistakes])
        print("="*50)
        print("💀 YOU LOST! 💀")
        print(f"The word was: {word_to_guess}")
        print("="*50)

def main():
    """Main function to run the game."""
    while True:
        play_game()
        
        # Ask to play again
        play_again = input("\nDo you want to play again? (yes/no): ").strip().lower()
        if play_again not in ["yes", "y"]:
            print("\nThanks for playing! Goodbye! 👋")
            break

if __name__ == "__main__":
    main()
    
