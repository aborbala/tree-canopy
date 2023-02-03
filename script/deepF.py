import logging

logging.basicConfig(filename='C:/MyLog.log', level=logging.DEBUG)

try:
    import deepforest as dp
    model = dp.main.deepforest()
    sample_image = dp.get_data("OSBS_029.png")
    print('YEEEEEEEEEEEEEEEEY')
    print(sample_image)

except Exception as e:
    logging.info(e)