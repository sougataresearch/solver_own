import sys
sys.path.append(r'C:\Program Files\Lumerical\v241\api\python')
import lumapi
from collections import OrderedDict
import math


class Substrate():
    def __init__(self) -> None:
        # 可以存一些預設的參數
        self.material = "Si (Silicon) - Palik" # 要非常精確，用GUI 左下角輸入 ? set("meterial") 用複製的過來
        self.x = 0
        self.y = 0
        self.z = 0
        self.x_span = 30e-6
        self.y_span = 20e-6
        self.y_min = -500e-6 # um 
        self.z_span = 30e-6
        
    
    def generate(self, app_handle: lumapi.FDTD) -> None: #這邊的: 是指定資料類型 e.g. float 
        '''
        建立substrate
        '''
        
        substrate_prop = OrderedDict([
                                    ("name", "substrate"),
                                    ("material", self.material),
                                    ("x", self.x), # 單位meter
                                    ("y", self.y), # 單位meter
                                    ("z", self.z), # 單位meter
                                    ("x span", self.x_span), # 單位meter 
                                    ("y span", self.y_span), # 單位meter
                                    ("y min", self.y_min), # 單位meter
                                    ("z span", self.z_span) # 單位meter 
                                    ])    
    
        app_handle.addrect(properties = substrate_prop)

    def get_surface_y(self):
        return self.y_span/2
 
class Grating():
    def __init__(self, substrate_Z_span, surface_y_position, x = 0, y = 0, z = 0) -> None:
        # 可以存一些預設的參數
        self.x = x
        self.y = y
        self.z = z
        
        
        
        self.material = "Si (Silicon) - Palik" # 要非常精確，用GUI 左下角輸入 ? set("meterial") 用複製的過來
        self.grating_number = 5
        self.ttcd = 1.5 * (10**-6) # meter
        self.tcd = 2 * (10**-6) # meter
        # self.bbcd = 3 * (10**-6) # meter
        self.bcd = 2.5 * (10**-6) # meter
        self.depth = 4 * (10**-6) # meter
        self.depth_of_top = 0.3 * (10**-6) # meter
        # self.depth_of_base = 0.5 * (10**-6) # meter
        self.spacing = 2 * (10**-6) # meter


        self.z_span = substrate_Z_span
        self.y_compensation = surface_y_position
        self.short_axis = 0.2 * (10**-6) # meter         
        
        
        
    def generate(self, app_handle: lumapi.FDTD) -> None: #這邊的: 是指定資料類型 e.g. float 
        app_handle.addstructuregroup() # 建立結構群組
        app_handle.set("name","grating")
        app_handle.set("x", self.x)
        app_handle.set("y", self.y)
        app_handle.set("z", self.z)

        # 在 group structure 裡面新增變數 
        '''
        
        adduserprop("property name", type, value);
        
        ref: https://optics.ansys.com/hc/en-us/articles/360034928733-adduserprop-Script-command
        Adds a user property to a selected structure group. The name is set to "property name". The type is an integer from 0 to 6. The corresponding variable types are

        0 - number

        1 - Text

        2 - Length

        3 - Time

        4 - Frequency

        5 - Material

        6 - Matrix

        The value of the new user property is set to value.
        '''
 
        app_handle.adduserprop("grating_number", 0, self.grating_number) 
        app_handle.adduserprop("material", 5, self.material)

        app_handle.adduserprop("ttcd", 2, self.ttcd)
        app_handle.adduserprop("tcd", 2, self.tcd)
        # app_handle.adduserprop("bbcd", 2, self.bbcd)
        app_handle.adduserprop("bcd", 2, self.bcd)
        app_handle.adduserprop("depth", 2, self.depth)
        app_handle.adduserprop("DoT", 2, self.depth_of_top)
        # app_handle.adduserprop("DoB", 2, self.depth_of_base) 
        app_handle.adduserprop("spacing", 2, self.spacing)

        app_handle.adduserprop("zSpan", 2, self.z_span)
        app_handle.adduserprop("yCompensation", 2, self.y_compensation)
        app_handle.adduserprop("shortAxis", 2, self.short_axis)
        
        
        ## 這個script 裡面的內容建議在 GUI 裡面測完再貼過來，不然很難 De-bug
        app_handle.set("script","""
                       
                        deleteall;

                        # --- 參數設定 ---
                        period = tcd + spacing;

                        # --- 頂點定義 ---

                        etch_vertices = [-(bcd/2), yCompensation;
                                        bcd/2, yCompensation; 
                                        tcd/2, (depth) + yCompensation; 
                                        -(tcd/2), (depth) + yCompensation];

                        upper_vertices = [-(tcd/2), (depth) + yCompensation;
                                        tcd/2, (depth) + yCompensation;
                                        ttcd/2, (depth + DoT) + yCompensation;
                                        -(ttcd/2), (depth + DoT) + yCompensation];

                        # --- 1. 建立中心物件 (Index 0) ---
                        # Base 0
                        # addpoly; set("name", "base_0");
                        # set("vertices", base_vertices);
                        # set("z span", zSpan); set("material", material);

                        # Etch 0
                        addpoly; set("name", "etch_0");
                        set("vertices", etch_vertices);
                        set("z span", zSpan); set("material", material);

                        # Upper 0
                        addpoly; set("name", "upper_0");
                        set("vertices", upper_vertices);
                        set("z span", zSpan); set("material", material);

                        # --- 2. 使用單一迴圈進行整組複製 ---
                        for(i=1; i<(grating_number+3)/2; i=i+1){
                            # 向右複製 (+x方向)
                            # unselectall;
                            # select("base_0");  copy(period*i, 0, 0); set("name", "base_" + num2str(i));
                            unselectall;
                            select("etch_0");  copy(period*i, 0, 0); set("name", "etch_" + num2str(i));
                            unselectall;
                            select("upper_0"); copy(period*i, 0, 0); set("name", "upper_" + num2str(i));

                            # 向左複製 (-x方向)
                            # unselectall;
                            # select("base_0");  copy(-period*i, 0, 0); set("name", "base_" + num2str(-i));
                            unselectall;
                            select("etch_0");  copy(-period*i, 0, 0); set("name", "etch_" + num2str(-i));
                            unselectall;
                            select("upper_0"); copy(-period*i, 0, 0); set("name", "upper_" + num2str(-i));
                        }
                       
                       
                       
                       
                       """)
        
    def get_surface_y(self):
        return self.y_compensation + self.depth + self.depth_of_top # 基材頂部 + grating 高度 + 圓頂短半軸
    
    def get_period(self):
        return self.tcd + self.spacing

class RCWARegion():
    def __init__(self, x_span = 10e-6, z_span = 10e-6, surface_y_position= 10e-6, interface = [("::model::substrate", "max", 1)]):
        
        #interface = [("::model::substrate", "max", 1)] 這是選定substrate 然後切一層
        
        #設定RCWA region 的中心以及大小
        self.x = 0
        self.z = 0 
        self.x_span = x_span # 單一週期的尺寸
        self.z_span = z_span # 也是單一週期的尺寸
        
        # RCWA 的Y 最高跟最低
        self.y_max = surface_y_position + 5e-6
        self.y_min = surface_y_position - 550e-6

        # RCWA 的rectangular K-vector 設定
        self.ku_vector = 15
        self.kv_vector = 1  
        # K-vector = 2N + 1
        self.k_vector = 101
        
        '''
        ref: https://optics.ansys.com/hc/en-us/articles/12959229278611-RCWA-Solver-Simulation-Object
        Solver Tab
        
        懶人包: 雙軸的K-vector 相乘，但不一定兩邊是一樣階數，總階數會很接近
        '''
        self.interface = interface 
        self.minimum_wavelength = 0.20e-6 #最大波長 250 nm
        self.maximum_wavelength = 0.40e-6 #最小波長 450 nm
        self.angle_theta = 0 # 入射角
        self.angle_phi = 0 # 方位角
        self.frequency_points = 100 # 波長要採樣多少點
        self.use_wavelength_spacing = False
        
    def generate(self, app_handle:lumapi.FDTD):

        solver_region_prop = OrderedDict([("propagation axis", "y"),
                                          ("propagation direction", "backward"),
                                          
                                          ("x", self.x),
                                          ("z", self.z),
                                          ("x span", self.x_span),
                                          ("z span", self.z_span),
                                          ("y max", self.y_max),
                                          ("y min", self.y_min),
                                          
                                          # rectangular k-vector settings
                                          ("k vectors domain", "rectangular"),
                                          ("max number ku", self.ku_vector),
                                          ("max number kv", self.kv_vector),
                                          # circular k-vector settings
                                        #   ("k vectors domain", "circular"),
                                        #   ("max number k vectors", self.k_vector),
                                          
                                          ("maximum wavelength", self.maximum_wavelength),
                                          ("minimum wavelength", self.minimum_wavelength),
                                          
                                          ("interface position", "reference"),  
                                          
                                          # reference: 基於結構的介面的 max or min
                                          # absolute: 基於座標下去切 
                                          
                                          ("interface reference positions", self.interface), #RCWA 要怎麼切層數
                                          ("report grating orders", True), # 看要求得什麼結果
                                          ("report grating characterization", True),
                                          # 設定confomal 1 
                                          ("mesh refinement", "conformal variant 1"),
                                          
                                          #波長採樣的方法
                                          ("use wavelength spacing", self.use_wavelength_spacing),
                                          ("frequency points", self.frequency_points),
                                          
                                          # 設定入射角
                                          ("angle theta", self.angle_theta),
                                          ("angle phi", self.angle_phi)])
        
        '''
        ref: https://optics.ansys.com/hc/en-us/articles/12959229278611-RCWA-Solver-Simulation-Object
                                          看更多: In ANALYSIS mode 可以分析的東西
        '''
        
        # 新增一個RCWA 的物件
        app_handle.addrcwa()
        
        # 逐個填進去
        for i in solver_region_prop:
            app_handle.set(i, solver_region_prop[i])
        
if __name__ == '__main__':
    ## 開FDTD的程式
    #fdtd_handle = lumapi.FDTD(filename=r'./Grating sample_t.fsp', hide= False) # hide 是用來決定要不要顯示GPI 的
    fdtd_handle = lumapi.FDTD(hide= False) # hide 是用來決定要不要顯示GPI 的
    
    # 創物件
    substrate = Substrate()
    substrate.generate(fdtd_handle)
    
    # 創grating 
    grating = Grating(substrate.z_span, substrate.get_surface_y(), substrate.x, substrate.y, substrate.z)
    grating.generate(fdtd_handle)
    
    #創RCWA
    solver_region_x_span = grating.get_period()
    
    ## 這邊的名子要回去FDTD GUI對造一下
    # rcwa_interface = [["::model::grating", "max", 4],
    #                   ["::model::grating::oval(0)", "min", 1],
    #                   ["::model::substrate", "max", 1],
    #                   ["::model::substrate", "min", 1]]
    
    
    rcwa_interface = [["::model::grating","max", 3],
                      ["::model::grating::etch_0", "max", 20],
                    #   ["::model::grating::base_0", "max", 5],
                      ["::model::substrate", "max", 1],
                      ["::model::substrate", "min", 1]]
    solver_region  = RCWARegion(x_span = solver_region_x_span, # 一個週期尺寸
                                z_span= solver_region_x_span, # 也是一個周期的尺寸
                                surface_y_position = grating.get_surface_y(), # 取得結構的頂部位置
                                interface = rcwa_interface # 切幾層
                                )
    
    solver_region.generate(fdtd_handle)
    
    fdtd_handle.save("grating_sim_5CD.fsp")
    

    
    
    _ = input('按任何鍵關閉')
    fdtd_handle.close() # 正常關閉 FDTD 程式
    print('關閉FDTD')