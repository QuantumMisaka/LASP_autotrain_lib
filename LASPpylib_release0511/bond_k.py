class bond(object):
    def __init__(self, ele1,ele2,dis):
        self.ele1 = ele1
        self.ele2 = ele2
        self.length = dis

    def judge_bondorder(self):
        ele1, ele2 = sorted([self.ele1, self.ele2])
        # for C H O only
        if ele1 == 1:
            if   ele2 == 1:
                if self.length < 0.84: return 1
                else:                  return 0
            elif ele2 == 5:
                if self.length < 1.35: return 1
                else:                  return 0
            elif ele2 == 6:
                if self.length < 1.19: return 1
                else:                  return 0
            elif ele2 == 7:
                if self.length < 1.11: return 1
                else:                  return 0
            elif ele2 == 8:
                if self.length < 1.06: return 1
                else:                  return 0
            elif ele2 == 9:
                if self.length < 1.02: return 1
                else:                  return 0
            elif ele2 == 14:
                if self.length < 1.58: return 1
                else:                  return 0
            elif ele2 == 15:
                if self.length < 1.54: return 1
                else:                  return 0
            elif ele2 == 16:
                if self.length < 1.44: return 1
                else:                  return 0
            elif ele2 == 17:
                if self.length < 1.37: return 1
                else:                  return 0
            elif ele2 == 32:
                if self.length < 1.63: return 1
                else:                  return 0
            elif ele2 == 34:
                if self.length < 1.56: return 1
                else:                  return 0
            elif ele2 == 35:
                if self.length < 1.51: return 1
                else:                  return 0
            elif ele2 == 50:
                if self.length < 1.8:  return 1
                else:                  return 0
            elif ele2 == 52:
                if self.length < 1.8:  return 1
                else:                  return 0
            elif ele2 == 53:
                if self.length < 1.71: return 1
                else:                  return 0

        # boron
        elif ele1 == 5 :
            if ele2 == 6:
                if self.length < 1.66: return 1
                else:                  return 0
            elif ele2 == 17:
                if self.length < 1.85: return 1
                else:                  return 0
        # carbon
        elif ele1 == 6 :
            if ele2 == 6:
                if   self.length < 1.25: return 3
                elif self.length < 1.44: return 2
                elif self.length < 1.65: return 1
                else:                    return 0
            elif ele2 == 7:
                if   self.length < 1.26: return 3
                elif self.length < 1.40: return 2
                elif self.length < 1.58: return 1
                else:                    return 0
            elif ele2 == 8:
                if   self.length < 1.15: return 3
                elif self.length < 1.3 : return 2
                elif self.length < 1.55: return 1
                else:                    return 0
            elif ele2 == 9:
                if self.length < 1.48: return 1
                #if self.length < 1.52: return 1
                else:                  return 0
            elif ele2 == 14:
                if self.length < 1.96: return 1
                else:                  return 0
            elif ele2 == 15:
                if self.length < 1.97: return 1
                else:                  return 0
            elif ele2 == 16:
                if   self.length < 1.7 : return 2
                elif self.length < 1.92: return 1
                else:                    return 0
            elif ele2 == 17:
                if self.length < 2.09: return 1
                else:                  return 0
            elif ele2 == 32:
                if self.length < 2.05: return 1
                else:                  return 0
            elif ele2 == 35:
                if self.length < 2.04: return 1
                else:                  return 0
            elif ele2 == 53:
                if self.length < 2.24: return 1
                else:                  return 0
        # nitrogen
        elif ele1 == 7 :
            if ele2 == 7:
                if   self.length < 1.2:  return 3
                elif self.length < 1.35: return 2
                elif self.length < 1.55: return 1
                else:                    return 0
            elif ele2 == 8:
                if   self.length < 1.24: return 2
                elif self.length < 1.56: return 1
                else:                    return 0
        # oxygen
        elif ele1 == 8 :
            if ele2 == 8:
                if   self.length < 1.3 : return 2
                elif self.length < 1.58: return 1
                else:                    return 0
            elif ele2 == 14:
                if self.length < 1.73: return 1
                else:                  return 0
            elif ele2 == 15:
                if   self.length < 1.48: return 2
                elif self.length < 1.73: return 1
                else:                    return 0
            elif ele2 == 16:
                if self.length < 1.53: return 1
                else:                  return 0

        return 0
