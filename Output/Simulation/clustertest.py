
from cc3d import CompuCellSetup
        

from WallTest2Steppables import WallTest2Steppable

CompuCellSetup.register_steppable(steppable=WallTest2Steppable(frequency=1))


CompuCellSetup.run()
