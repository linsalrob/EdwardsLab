"""
An answer to http://fivethirtyeight.com/features/how-long-will-your-smartphone-distract-you-from-family-dinner/
"""

from random import choice

if __name__ == '__main__':
    tasks = [1,2,3,4,5]
    total = []
    for i in range(1000):
        mine = choice(tasks)
        sisters = choice(tasks)
        while mine != sisters:
            mine += choice(tasks)
            sisters += choice(tasks)
        total.append(mine)
    print(1.0 * sum(total)/len(total))

